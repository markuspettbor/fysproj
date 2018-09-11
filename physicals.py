import numpy as np
import numtools as nt
import prob_dist as pd
import operator

class Gas:
    def __init__(self, num_particles, temperature, low = 0, high = 1, mass = 0, radius  = 0):
        self.num_particles = num_particles
        self.temperature = temperature
        self.position = np.array([])
        self.velocity = np.array([])
        self.createParticles(num_particles, mass, low, high, radius)

    def createParticles(self, new, mass, low, high, radius = 0):
        boltzmax = pd.BoltzmannMaxwell(self.temperature, new)
        self.position = np.random.uniform(low, high, size = (3, new))
        self.velocity = boltzmax.distribution(size = (3, new))

    def addParticles(self):
        # Syntax: np.append(vector)
        pass


class Wall():
    def __init__(self, normal_vector, axis, sign, center, hole_width = 0, corners = None, molecule_moment = 0):
        self.normal_vector = normal_vector
        self.axis = axis
        self.sign = sign
        self.center = center
        self.corners = corners
        self.hole_width = hole_width
        self.unit_normal = normal_vector/np.linalg.norm(normal_vector)
        self.escaped_particles = 0
        self.escaped_velocity = 0

    def check_collision(self, position, velocity):
        velocity[self.axis] = np.where(self.boundary(position, velocity),\
                            -velocity[self.axis], velocity[self.axis])
        return velocity

    def boundary(self, position, velocity):
        '''
        Assumes position is a numpy array containing positions (x,y,z),
        with dimension (3, N).
        First checks if any particle is outside the wall along a given axis.
        If there is a hole, and particles are outside, it checks to see if those
        particles are within the bounds of the hole.
        '''
        outside = self.sign*position[self.axis] > self.sign*self.center[self.axis]
        # Find out which particles are outside. Then, check among those that are outside, if they are in hole!
        match = np.any(outside)
        if self.hole_width != 0 and match:
            mask_index = np.arange(3) != self.axis # retrieve relevant axes
            grid = np.where(outside, position[mask_index], False)
            in_hole = np.abs(self.center[mask_index] - grid.transpose()) <= self.hole_width/2
            in_hole = in_hole.transpose()
            esc = in_hole[0]*in_hole[1]*outside
            self.escaped_particles += np.count_nonzero(esc)
            self.escaped_velocity += (np.sum(np.abs(velocity[self.axis, esc]))) #norm is not work
        return outside

    def reset_escaped_particles(self):
        self.escaped_particles = 0
        self.escaped_velocity = 0


class Box:
    def __init__(self, dx, dy, dz, width, hole_index = 5, hole_width = 200):
        '''
        Creates the simplest of all boxes.
        '''
        normals = np.array([[1,0,0], [-1, 0, 0], [0, 1, 0],\
                            [0, -1, 0], [0, 0, 1], [0, 0, -1]])
        axes = np.array([0, 0, 1, 1, 2, 2])
        signs = np.array([1, -1, 1, -1, 1, -1])
        x1 = np.array([dx + width/2, dy, dz])
        x2 = np.array([dx - width/2, dy, dz])
        x3 = np.array([dx, dy + width/2, dz])
        x4 = np.array([dx, dy - width/2, dz])
        x5 = np.array([dx, dy, dz + width/2])
        x6 = np.array([dx, dy, dz - width/2])
        centers = [x1, x2, x3, x4, x5, x6]
        self.walls = []
        holes = np.zeros(6)
        holes[hole_index] = hole_width
        self.make_box(normals, axes, signs, centers, holes)

    def add_wall(self, normal, axis, sign, center, hole):
        return Wall(normal, axis, sign, center, hole_width = hole)

    def make_box(self, norms, axes, signs, centers, holes):
        for n, ax, sign, cntr, hole in zip(norms, axes, signs, centers, holes):
            self.walls.append(self.add_wall(n, ax, sign, cntr, hole))


if __name__ == '__main__':
    print('main')

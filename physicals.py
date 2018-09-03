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
    def __init__(self, normal_vector, axis, sign, center, corners = None, hole_width = 0, molecule_moment = 0):
        self.normal_vector = normal_vector
        self.axis = axis
        self.sign = sign
        self.center = center
        self.corners = corners
        self.hole_width = hole_width
        self.unit_normal = normal_vector/np.linalg.norm(normal_vector)
        self.escaped = 0
        self.dp = 0
        self.p = molecule_moment

    def check_collision(self, position, velocity):
        velocity[self.axis] = np.where(self.boundary(position),\
                            -velocity[self.axis], velocity[self.axis])
        return velocity

    def boundary(self, position):
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
            self.escaped += np.count_nonzero(esc)
            print(self.escaped)
        return outside


if __name__ == '__main__':
    print('main')

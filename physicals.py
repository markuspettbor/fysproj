import numpy as np

class Gas:
    def __init__(self, num_particles, temperature, positions, radius  = 0):
        self.num_particles = num_particles
        self.temperature = temperature
        self.positions = positions

    def addParticles(self, new_particles, mass, radius = 0):
        pass

class Wall:
    def __init__(self, w, h, normal_vector, position):
        '''
        w, h, is the width and height of a rectangular wall.
        Assumes normal_vector is a numpy array on the form [x, y, z]
        The normal_vector determines the orientation of the wall relative
        to the center of a box.
        '''
        self.w = w
        self.h = h
        self.normal_vector = normal_vector
        self.position = position

    def __call__(self):
        print('wall coordinates returned')

    def get_corners(self):
        '''
        Finds corners in a rectangle of width w, height h, for a given
        normal_vector. Idea: Find a vector in the plane, cross that vector
        with the normal vector, and get a vector perpendicular to the vector
        in the plane. The two can then be used to span a rectangle.

        For future: Should redo entire code. Should also add option
        to build a wall from an existing one.
        '''
        a, b, c = self.normal_vector
        x0, y0, z0 = self.position
        p = np.array([0, 0, a/c*x0 + b/c*y0 + z0])
        tryvec = p - self.position
        othervec = np.cross(tryvec, self.normal_vector)
        res1 = self.position + tryvec/np.linalg.norm(tryvec)*self.w + othervec/np.linalg.norm(othervec)*self.h
        res2 = self.position + tryvec/np.linalg.norm(tryvec)*self.w
        res3 = self.position + othervec/np.linalg.norm(othervec)*self.h
        return res1, res2, res3, self.position


def build_the_wall():
    wall = Wall(10, 5, np.array([0, 0, 1]), np.array([1,1,1]))
    wall.get_corners()

build_the_wall()

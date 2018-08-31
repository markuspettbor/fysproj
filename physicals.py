import numpy as np
import numtools as nt

class Gas:
    def __init__(self, num_particles, temperature, positions, radius  = 0):
        self.num_particles = num_particles
        self.temperature = temperature
        self.positions = positions

    def addParticles(self, new_particles, mass, radius = 0):
        pass

class Wall:
    def __init__(self, w, h, normal_vector, position, neighbours = None):
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
        p = np.zeros(3)
        try:
            p1 = b/a*y0 + c/a*z0 + x0
        except ValueError:
        p2 = a/b*x0 + c/b*z0 + y0
        p3 = a/c*x0 + b/c*y0 + z0

        for point, index in zip((p1, p2, p3), (0, 1, 2)):
            if point is not 0:
                p[index] = point
                break

        u1 = (p - self.position)nt.norm
        u1 = r1/nt.norm()
        r2 = np.cross(r1, self.normal_vector)
        c1 = self.position + r1/nt.norm(r1)*self.w/2 + r2/nt.norm(r2)*self.h/2
        c2 = self.position + r1/nt.norm(r1)*self.w/2 - r2/nt.norm(r2)*self.h/2
        c3 = self.position - r1/nt.norm(r1)*self.w/2 + r2/nt.norm(r2)*self.h/2
        c4 = self.position - r1/nt.norm(r1)*self.w/2 - r2/nt.norm(r2)*self.h/2
        return np.array([c1, c2, c3, c4])


def build_the_wall():
    wall = Wall(100, 100, np.array([0,0, 1]), np.array([1000, 1000, 0]))
    print(wall.get_corners())

build_the_wall()

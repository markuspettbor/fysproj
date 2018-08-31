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
        pass



def build_the_wall():
    pass

build_the_wall()

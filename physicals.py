import numpy as np
import numtools as nt
import prob_dist as pd

import matplotlib.pyplot as plt

class Gas:
    def __init__(self, num_particles, temperature, mass = 0, radius  = 0):
        self.num_particles = num_particles
        self.temperature = temperature
        self.position = np.zeros((1, 3))
        self.velocity = np.zeros((1, 3))
        self.addParticles(num_particles, mass, radius)

    def addParticles(self, new, mass, radius = 0):
        boltzmax = pd.BoltzmannMaxwell()
        np.append(self.position, np.random.uniform(size = (new, 3)))
        np.append(self.velocity, boltzmax.distribution(size = (new, 3) ) )
        print(self.velocity)

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


def test():
    h2 = Gas(10, 5000)
    print(h2.num_particles)
    plt.plot(h2.velocity)
    plt.show()

test()

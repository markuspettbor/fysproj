import numpy as np

class Gas:
    def __init__(self, num_particles, temperature, positions, radius  = 0):
        self.num_particles = num_particles
        self.temperature = temperature
        self.positions = positions

    def addParticles(self, new_particles, mass, radius = 0):
        pass

class Wall:
    def __init__(self, w, h, fixed):
        self.w = w
        self.h = h
        self.fixed = fixed

    def __call__(self):
        print('wall coordinates returned')

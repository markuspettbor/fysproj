import numpy as np

class Gas:
    def __init__(self, num_particles, temperature, positions):
        self.num_particles = num_particles
        self.temperature = temperature
        self.positions = positions

    def addParticles(self, new_particles, mass):
        pass

class Wall:
    def __init__(self):
        pass

    def __call__(self):
        print('wall coordinates returned')

class Engine:
    def __init__(self, length, width, height):
        pass

    def make_wall(self):
        pass

    def detect_collision(self):
        pass

#[x, y, vx, vy, mass, radius, force]

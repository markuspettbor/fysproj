import numpy as np

class SolarSystem:
    def __init__(self, num_planets, num_suns = 1):
        self.num_planets = num_planets
        self.planets = {}
        self.suns = {}

    def addPlanet(self, planet, name):
        self.planets[name] = planet

    def addSun(self, sun, name):
        self.suns[name] = sun

    def gravity(self, m1, m2, r):
        pass

class Planet:
    def __init__(self, mass, radius, initial_position):
        self.mass = mass
        self.radius = radius
        self.initial_position = initial_position

class Sun:
    def __init__(self, mass, radius, initial_position, temperature):
        self.mass = mass
        self.radius = radius
        self.initial_position = initial_position
        self.temperature = temperature

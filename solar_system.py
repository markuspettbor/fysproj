import numpy as np

class SolarSystem:
    def __init__(self, num_planets, num_suns = 1):
        self.num_planets = num_planets
        self.planets = []
        self.suns = []

    def addPlanet(self, planet):
        self.planets.append(planet)

    def addSun(self, sun):
        self.suns.append(sun)

    def gravity(self, m1, m2, r):
        pass

class Planet:
    def __init__(self, mass, radius, position, velocity, name):
        self.mass = mass
        self.radius = radius
        self.position = position
        self.velocity = velocity
        self.name = name

class Sun:
    def __init__(self, mass, radius, position, temperature, name):
        self.mass = mass
        self.radius = radius
        self.position = position
        self.temperature = temperature
        self.name = name

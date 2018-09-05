import numpy as np
import variables

class Chamber:
    def __init__(self):
        pass

class Engine:
    def __init__(self, fuel_mass, mass):#, rocket):
        '''
        rocket is an object of the rocket to be simulated
        needs this for its total mass
        '''
        #self.rocket = rocket
        self.fuel_mass = fuel_mass #mass of fuel   given by rocket?
        self.mass = mass #mass of engine           given by rocket?
        self.molecule_mass = variables.molecule_mass

    def particles_escaped(self, particles, velocity, delta_time):
        self.fuel_mass += particles*variables.molecule_mass*(-1)
        delta_momentum = variables.molecule_mass*velocity
        force = delta_momentum/delta_time
        return force

    def boost(self, rocket_thrust_force, fuel_consumption, initial_rocket_mass, delta_speed):
        while speed > delta_speed:
            pass
    '''
    can either return force, or save force as self.force and call
    engine to get the force
    '''

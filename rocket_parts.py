import numpy as np
import variables
import numtools as nt

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
        fuel_used = particles*variables.molecule_mass
        delta_momentum = variables.molecule_mass*velocity
        force = delta_momentum/delta_time
        return force, fuel_used


    '''
    can either return force, or save force as self.force and call
    engine to get the force
    '''

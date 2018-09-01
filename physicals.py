import numpy as np
import numtools as nt
import prob_dist as pd
<<<<<<< HEAD
import matplotlib.pyplot as plt
import operator

=======
>>>>>>> ebc446e063dbb53a42014bf461ce8caa33987c02

class Gas:
    def __init__(self, num_particles, temperature, mass = 0, radius  = 0):
        self.num_particles = num_particles
        self.temperature = temperature
        self.position = np.array([])
        self.velocity = np.array([])
        self.createParticles(num_particles, mass, radius)

    def createParticles(self, new, mass, radius = 0):
        boltzmax = pd.BoltzmannMaxwell(self.temperature, new)
        self.position = np.random.uniform(size = (new, 3))
        self.velocity = boltzmax.distribution(size = (new, 3))

    def addParticles(self):
        # Syntax: np.append(vector)
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


    def get_wall_normal(self, corners, origin, index = 1): #name shoudl be initialise wall??
        '''
        corners given as
        index is either -1 or 1, index sier om normalvektor skal inn eller ut av origo
        NOTAT: om vinkelen mellom normalvektoren og en possisjonsvektor på planet er
        mindre enn 90 grader må den snus (*-1), gitt at dette punktet ikke ligger i origo
        '''
        try: #not sure if non-self is much faster
            corners = np.array(corners)
            origin = np.array(origin)
            index = int(index)
            self.corners = corners
            self.origin = origin
        except ValueError:
            print('corners are not one array, or index is not int')

        p0p1 = corners[1] - corners[0]#p1 - p0, vektor fra hjørne 0 til hjørne 1
        p0p2 = corners[2] - corners[0]#p2 - p0, vektor fra hjørne 0 til hjørne 2
        n_temp = np.cross(p0p1, p0p2)
        p0_temp= corners[0] - origin #vector from origin to a point in the plane
        angle = angle_between(n_temp, p0_temp)
        if angle > np.pi/2:
            n = unit_vector(n_temp)*index #index bestemmer hva som er inn og ut av legemet
        elif angle < np.pi/2:
            n = unit_vector(n_temp)*index*(-1) #(-1) snur normalvektoren inn i legemet
        else:
            raise ValueError('the plane goes through the origin, not ok') #may want to do something hardcoded here
        return n

        #finn hvilken retning (x, y eller z) som normalvektoren øker fortest i
        #(største komponent), sett de 2 andre konstant og finn ut når den krysser
        #denne 'komponenten'
        index_n, value = max(enumerate(np.abs(n)), key=operator.itemgetter(1))
        #https://stackoverflow.com/questions/6193498/pythonic-way-to-find-maximum-value-and-its-index-in-a-list
        self.teller = np.sum(n*(corners[0]-origin) #ax0+by0+cz0
        self.index_n = index_n

    def gain_momentum(self, position, hole_hw):
        n = self.n
        i = 0
        for pos in position: #burde kunne parallelliseres litt
            if pos[2]) < -hw: #ser om den er utenfor bredden
                if abs(pos[1]<hole_hw) or abs(pos[0]<hole_hw):
                    return True
                else:


    def new_vel(self, position, velocity, hw): #hw is halfwidth, sym box, en vegg er 2 vegger om det gir mening, speilet om origo
        n = self.n
        index = self.index_n
        i = 0
        for pos, vel in zip(position, velocity): #burde kunne parallelliseres litt
            if abs(pos[index]) > hw: #ser om den er utenfor bredden
                if np.sign(pos[index]) != np.sign(n[index]):
                    velocity[i] = vel*n
                else:
                    velocity[i] = -vel*n
        i = i + 1
'''
    def get_new_velocity(self, position, velocity): #wall equation to be called
        index = self.index_n
        for pos, vel in zip(position, velocity): #burde kunne parallelliseres litt
        try:

            if index == 0: #normalvektor peker mest i x retning
                if
                index_is_zero(pos, vel)
            elif index == 1: #normalvektor peker mest i y retning
                index_is_one(pos, vel)
            elif index == 2: #normalvektor peker mest i z retning
                index_is_two(pos, vel)
            else:
                raise ValueError('n vector had no max value??')

    def index_is_zero(self, position, velocity): #then use this as function
        n = self.n
        pos_eq = self.teller/(n[1]*pos[1]+n[2]*pos[2])/n[0]
        if np.sign() =

    def index_is_one(self, position, velocity): #then use this as function
        n = self.n
        pos_eq = self.teller/(n[0]*pos[0]+n[2]*pos[2])/n[1]

    def index_is_two(self, position, velocity): #then use this as function
        n = self.n
        pos_eq = self.teller/(n[0]*pos[0]+n[1]*pos[1])/n[2]
'''



def build_the_wall():
    pass

build_the_wall()

def unit_vector(vector):
    """ https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    Returns the unit vector of the vector.  """
    # Please use norm or unit_vector functions from numtools module for beautifullnes
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    Returns the angle in radians between vectors 'v1' and 'v2':: """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

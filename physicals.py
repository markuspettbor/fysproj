import numpy as np
import numtools as nt
import prob_dist as pd

class Gas:
    def __init__(self, num_particles, temperature, mass = 0, radius  = 0):
        self.num_particles = num_particles
        self.temperature = temperature
        self.position = np.array([])
        self.velocity = np.array([])
        self.createParticles(num_particles, mass, radius)

    def createParticles(self, new, mass, radius = 0):
        boltzmax = pd.BoltzmannMaxwell()
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


    def get_wall_normal(self, corners, origin, index = 1):
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
            self.n = unit_vector(n_temp)*index #index bestemmer hva som er inn og ut av legemet
        elif angle < np.pi/2:
            self.n = unit_vector(n_temp)*index*(-1) #(-1) snur normalvektoren inn i legemet
        else:
            raise ValueError('the plane goes through the origin, not ok') #may want to do something hardcoded here
        return n

        def get_wall_equation(self, corners, origin): #might not want this totally seperate from get_wall_normal
            corners = self.corners
            origin = self.origin
            n = self.n


def build_the_wall():
    pass

build_the_wall()

def unit_vector(vector):
    """ https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    Returns the angle in radians between vectors 'v1' and 'v2':: """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

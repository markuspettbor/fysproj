import numpy as np
import numtools as nt
import prob_dist as pd
import operator

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
    def __init__(self, d, normal_axis, hole_width = None, neighbours = None):
        '''
        d is distance from origin, negative number for negative position along the given axis
        Assumes normal_axis is 0 for x axis, 1 for y axis and 2 for z axis
        The normal_vector determines the orientation of the wall relative
        to the center of a box.
        '''
        normal_vector = [0,0,0]
        normal_vector[normal_axis] = 1*np.sign(d)
        self.d = d
        self.normal_vector = normal_vector
        self.sign = np.sign(normal_vector[normal_axis])
        self.index = normal_axis
        try:
            make_hole(hole_width)
        except: #vet ikke om jeg kan droppe det etter try helt
            pass

    def get_corners(self):
        pass

    def make_hole(self, w):
        self.w = w #half width of hole

    def reset_escape(self):
        self.escape_n = 0
        #self.escape_vel = 0
    def boundary(self, x):
        out = self.sign*x[self.index] > self.sign*self.d
        try:
            if out:
                w = self.w
                '''
                kjører dette if-kaoset bare hvis en partikkel krasjer med en vegg med hul
                ingen lur løsning på dette dukket opp i hodet mitt
                er kanskje ikke noe mer effektivt enn å sjekke exit med en "detect_exit()"
                får heller ikke ut hvilken hastighet den forlater med her, da måtte boundary
                tatt imot en hastighet for hver partikken den sjekker, men en kan vel bruke
                statistikk for å estimere utgangshastighet, men vil ikke det :(
                reset_escape() må vel calles etter en vegg er ferdig sjekket. altså før neste tidsintervall begynner.
                '''
                if self.index == 0:
                    if abs(x[1]) > w or abs(x[2]) > w:
                        self.escape_n = self.escape_n + 1
                elif self.index == 1:
                    if abs(x[0]) > w or abs(x[2]) > w:
                        self.escape_n = self.escape_n + 1
                elif self.index == 2:
                    if abs(x[0]) > w or abs(x[1]) > w:
                        self.escape_n = self.escape_n + 1
        finally:
            return out

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
            self.n = nt.unit_vector(n_temp)*index #index bestemmer hva som er inn og ut av legemet
        elif angle < np.pi/2:
            self.n = nt.unit_vector(n_temp)*index*(-1) #(-1) snur normalvektoren inn i legemet
        else:
            raise ValueError('the plane goes through the origin, not ok') #may want to do something hardcoded here
        return self.n

        def get_wall_equation(self, corners, origin): #might not want this totally seperate from get_wall_normal
            corners = self.corners
            origin = self.origin
            n = self.n

def build_the_wall():
    x = Wall(-1, 1)
    b = x.boundary([2,2,0])
    print(b)
    #x.check_collision(1, 2)
if __name__ == '__main__':
    build_the_wall()
'''
def unit_vector(vector):
    """ https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    Returns the unit vector of the vector.  """
    # Please use norm or unit_vector functions from numtools module for beautifullnes
    return vector / np.linalg.norm(vector)
'''
def angle_between(v1, v2):
    """ https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    Returns the angle in radians between vectors 'v1' and 'v2':: """
    v1_u = nt.unit_vector(v1)
    v2_u = nt.unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

from PIL import Image
import numpy as np
import numtools as nt
import variables as vars
from numba import jit
from numpy.linalg import inv as numpy_inv
from scipy.interpolate import interp1d as numpy_interpolate

class StereographicProjection:
    def __init__(self, fov_phi, fov_theta, phi0, theta0, img = None):
        self.fov_phi = fov_phi
        self.fov_theta = fov_theta
        self.phi0 = phi0
        self.theta0 = theta0
        self.img = img

    def xmaxmin(self):
        return 2*np.sin(self.fov_phi/2)/(1 + np.cos(self.fov_phi/2))

    def ymaxmin(self):
        return -2*np.sin(self.fov_theta/2)/(1 + np.cos(self.fov_theta/2))

    def find_phi(self, x, y):
        rho = np.sqrt(x**2 + y**2)
        beta = 2*np.arctan(rho/2)
        phi = self.phi0 + np.arctan(x*np.sin(beta)/\
            (rho*np.sin(self.theta0)*np.cos(beta)-y*np.cos(self.theta0)*np.sin(beta)))
        return phi

    def find_theta(self, x,y):
        rho = np.sqrt(x**2 + y**2)
        beta = 2*np.arctan(rho/2)
        theta = self.theta0 - np.arcsin(np.cos(beta)*np.cos(self.theta0)\
                            + y/rho*np.sin(beta)*np.sin(self.theta0))
        return theta

    def make_rgb(self, x, y, dims):
        make_img = np.zeros(dims)
        phi, theta = self.find_phi(x,y), self.find_theta(x,y)
        for i in range(dims[0]):
            for j in range(dims[1]):
                idx = vars.solar_system.ang2pix(theta[i,j],phi[i,j])
                rgb = sky[idx, 2:]
                make_img[i, j] = rgb
        return make_img

    def best_fit(self, ref, image):
        '''
        Tries to find best fit for projection centred around
        angle phi.
        Assumes ref is a reference array corresponding to a
        2 pi panorama of night sky. Assumes image is a
        part of that panorama, of same dimension as ref
        Returns best angle, found using least squares approach,
        in degrees :'(
        '''
        width = image.shape[1]
        rads_per_pixel = self.fov_phi/width
        ref = np.concatenate((ref, ref[:, :width]), axis = 1)
        best = np.sum(nt.norm(ref[:, :width] - image)**2)
        best_col = 0

        for col in range(ref.shape[1]- width):
            section = ref[:, col: col + width]
            distance = np.sum(nt.norm(section - image)**2)
            if distance < best:
                best = distance
                best_col = col
        return np.degrees((best_col + width/2)*rads_per_pixel)


def vel_rel_star(Dlam, lam0):
    vr = Dlam/lam0*vars.c #nm/nm*m/s -> [m/s]
    return vr[:,0] + vr[:,1]

def shift_ref(p, d, convert = 'from'):
    '''convert 'to' or 'from' the given coordinate system'''
    if convert == 'from':
        M = np.array([np.cos(p), np.sin(p)]).transpose()
        M_inv = numpy_inv(M)
        return M_inv * d
    elif convert == 'to':
        M = np.array([np.cos(p), np.sin(p)]).transpose()
        return  M * d
    else:
        print('convert has to be either to, or from. Default is from')
        pass

def velocity_from_stars(lam_measured, lam0 = 656.3): #phi [rad], lam [nm]
    '''calculates velocity in cartesian coordinates given angles and delta lambdas'''
    ref_stars = np.array(vars.ref_stars) #phi in DEGREES, lam in nanometers
    phi_skew = np.array(ref_stars[:,0])*np.pi/180 #phi in RADIANS
    lam_skew = np.array(ref_stars[:,1])
    lam_delta = lam_skew - lam_measured #difference between reference shift and measured shift

    lam = shift_ref(phi_skew, lam_delta, 'from')
    vel_cart = vel_rel_star(lam, lam0)# [m/s]
    #print('velocity of satelite with respect to sun in cartesian coordinates', vel_cart)
    return vel_cart/vars.AU_tall*vars.year

#TRILATERATION
#time of measurement: t0
#list of meadured distances: [p0, p1... pn, star]

def position_from_objects(current_time, distances, xx):
    #print('DIST', distances)
    d = np.zeros(len(distances))
    d[0] = distances[-1]
    d[1:] = distances[:-1]
    #print('d', d)
    #xx = np.load('saved/saved_orbits/launch_resolution/pos.npy')
    p = xx[:,:,current_time].transpose()
    #print('p', p)
    #d = np.random.random(9)     #distances
    #p = np.random.random([9,2]) #positions

    x = np.array([])
    y = np.array([])
    print(x)
    count = 0
    ind = np.linspace(0,len(d)-1,len(d), dtype = 'int')
    for i in ind:
        for j in ind[ind != i]:
            for k in ind[(ind != i) * (ind != j)]:
                #print(i,j,k)
                count += 1
                a2 = p[j,0] - p[i,0]    #a corresponds to x positions of planets
                a3 = p[k,0] - p[i,0]
                b2 = p[j,1] - p[i,1]    #b corresponds to y positions of planets
                b3 = p[k,1] - p[i,1]
                #c is a constant depandant on x and y positions in addition to distances
                c2 = d[i]**2 - d[j]**2 - p[i,0]**2 - p[i,1]**2 + p[j,0]**2 + p[j,1]**2
                c3 = d[i]**2 - d[k]**2 - p[i,0]**2 - p[i,1]**2 + p[k,0]**2 + p[k,1]**2

                y = np.append(y, 1/2*(a2*c3 - a3*c2) / (a2*b3 - a3*b2))
                x = np.append(x, (c3 - 2*y[i]*b3) / (2*a3))
    print(np.array([x,y]))
    print(count, 'COUNT')
    x_avg = np.average(x)
    y_avg = np.average(y)
    pos = np.array([x_avg, y_avg])
    return pos

'''
# Generate reference
phi0s = np.arange(np.ceil(360/fov_phi_deg))*fov_phi_deg + fov_phi_deg/2
ref = np.zeros(pixel_img.shape)
pixelsperdeg = pixel_img.shape[1]/fov_phi_deg
def dobaz():
    # 640 pixels = 70 degs,
    for phi0, i in zip(phi0s, range(len(phi0s))):
        phi0_rad = np.radians(phi0)
        projection.phi0 = phi0_rad
        x = np.linspace(-x_max, x_max, pixel_img.shape[1])
        y = np.linspace(-y_max, y_max, pixel_img.shape[0])
        xx, yy = np.meshgrid(x, y)
        next = projection.make_rgb(xx, yy, pixel_img.shape)
        if phi0 == phi0s[0]:
            ref = next
        else:
            ref = np.concatenate((ref, next), axis = 1)
    print(ref.shape)
    ref = ref[230:250, :int(pixelsperdeg*360)]
    print(ref.shape)
    np.save('saved/saved_params/reference_sky_ex.npy', ref)
    saver = Image.fromarray(ref.astype(np.uint8))
    saver.save('cool.png')
#dobaz()
'''
# Note that ref is a slice of the night sky, 20 pixels high.
def test():
    ref_img = Image.open('images/sample0000.png')
    pixel_img = np.array(ref_img)
    best = find_angle(pixel_img)
    print('Expected value: 0 degrees. Estimated value:', best)

def find_angle(picture):
    fov_phi_deg = 70
    fov_phi = 2*np.pi/360*fov_phi_deg
    fov_theta = fov_phi
    phi0 = 0
    theta0 = np.pi/2
    projection = StereographicProjection(fov_phi, fov_theta, phi0, theta0)

    x_max = projection.xmaxmin()
    y_max = projection.ymaxmin()

    sky = np.load('saved/saved_params/himmelkule.npy')
    sol_system = vars.solar_system

    ref = np.load('saved/saved_params/reference_sky_ex.npy')

    return projection.best_fit(ref, picture[230:250])
    # Assumes picture is 480x640 image

if __name__ == '__main__':
    test()

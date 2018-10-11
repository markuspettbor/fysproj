from PIL import Image
import numpy as np
import numtools as nt
import variables as vars
from numba import jit

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
        Tries to find best match for angle phi.
        Assumes ref is reference array of size (num_imgs*x, y, 3)
        Assumes image is comparison image, of size (x, y, 3)
        Returns index of best fit, which corresponds to best angular fit, in degs.
        '''
        smallest = 0
        num_imgs = len(ref)
        distance = np.zeros(num_imgs - 1)
        for img, i in zip(ref, range(1, num_imgs)):
            distance[i-1] = np.sum(nt.norm(img - image))
        return distance

    def experimental(self, ref, image):
        '''
        Tries to find best match for angle phi.
        Assumes ref is reference array of size(360/FOV*x, y, 3)
        Assumes image is comparison image, of size(x, y, 3)
        Assumes image is fully contained in at most two FOVs in reference
        '''
        distance = np.sum(nt.norm(ref - image))
        sort_dist = np.sort(distance)
        best = sort_dist[-2:]
        best_args = np.argwhere(distance == best)
        total = np.sum(best)
        weight = np.sum(best_args)/total



ref_img = Image.open('images/sample0000.png')
pixel_img = np.array(ref_img)

fov_phi = 2*np.pi/360*70
fov_theta = fov_phi
phi0 = 0
theta0 = np.pi/2
projection = StereographicProjection(fov_phi, fov_theta, phi0, theta0)

x_max = projection.xmaxmin()
y_max = projection.ymaxmin()

print('Max/min x: +- ', x_max, '\nMax/min y: +- ', y_max)
print('Full range x: ', 2*x_max, '\nFull range y: ', -2*y_max)
sky = np.load('saved/saved_params/himmelkule.npy')
sol_system = vars.solar_system

# Generate reference image of sky, hopefully only once.
'''
x_max = projection.xmaxmin()
y_max = projection.ymaxmin()
phi0s = np.linspace(0, 359, 360)*2*np.pi/360
ref = np.zeros((360, pixel_img.shape[0], pixel_img.shape[1], 3))
@jit
def dobaz():
    for phi0, angle in zip(phi0s, range(360)):
        print(angle)
        projection.phi0 = phi0
        x = np.linspace(-x_max, x_max, pixel_img.shape[1])
        y = np.linspace(-y_max, y_max, pixel_img.shape[0])
        xx, yy = np.meshgrid(x, y)
        ref[angle] = projection.make_rgb(xx, yy, pixel_img.shape)
    np.save('saved/saved_params/reference_sky3.npy', ref)
'''
reference = np.load('saved/saved_params/reference_sky3.npy')
def test():
    for i in range(10):
        minis = projection.best_fit(reference, reference[i])
        print(np.argmin(minis))

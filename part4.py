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
        best = np.sum(nt.norm(ref[:, :width] - image))
        best_col = 0

        for col in range(ref.shape[1]- width):
            section = ref[:, col: col + width]
            distance = np.sum(nt.norm(section - image))
            if distance < best:
                best = distance
                best_col = col
        return np.degrees((best_col + width/2)*rads_per_pixel)

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

# Generate reference for experimental

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

    fov_phi_deg = 70
    fov_phi = 2*np.pi/360*fov_phi_deg

    fov_theta = fov_phi
    phi0 = 0
    theta0 = np.pi/2
    projection = StereographicProjection(fov_phi, fov_theta, phi0, theta0)

    x_max = projection.xmaxmin()
    y_max = projection.ymaxmin()

    print('Max/min x: +- ', x_max, '\nMax/min y: +- ', y_max)
    print('Full range x: ', 2*x_max, '\nFull range y: ', -2*y_max)
    sky = np.load('saved/saved_params/himmelkule.npy')

    ref = np.load('saved/saved_params/reference_sky_ex.npy')
    width = 0 #int(ref.shape[1]/2-320)
    best = projection.best_fit(ref, ref[:, width:width + 640])
    print('Expected value: 35 degrees. Estimated value:', best)

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

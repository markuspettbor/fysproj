from PIL import Image
import numpy as np
import variables as vars

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
print('Full range x: ', 2*x_max, '\nFull range y: ', 2*y_max)
sky = np.load('saved/saved_params/himmelkule.npy')
sol_system = vars.solar_system

x = np.linspace(-x_max, x_max, pixel_img.shape[1])
y = np.linspace(-y_max, y_max, pixel_img.shape[0])

xx, yy = np.meshgrid(x, y)

phi, theta = projection.find_phi(xx,yy), projection.find_theta(xx,yy)

make_img = np.zeros(pixel_img.shape)
for i in range(480):
    for j in range(640):
        idx = sol_system.ang2pix(theta[i,j],phi[i,j])
        rgb = sky[idx, 2:]
        make_img[i, j] = rgb
saveimg = Image.fromarray(make_img.astype('uint8'))
saveimg.save('cool.png')

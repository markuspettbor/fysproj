from AST2000SolarSystem import AST2000SolarSystem
import matplotlib.pyplot as plt
import numpy as np
import numtools as nt
import prettyplot as pretty
import variables as vars

user = vars.seed
solar_system = vars.solar_system
n = vars.n
x0 = vars.x0
y0 = vars.y0
vx0 = vars.vx0
vy0 = vars.vy0
a = vars.a
e = vars.e
theta0 = vars.theta0
psi0 = vars.psi0
radius = vars.radius
m_star = vars.m_star
m = vars.m
G = vars.G

k = 1000
thetas = np.array([np.linspace(theta, theta + 2*np.pi, k) for theta in theta0])
thetas = np.transpose(thetas)
r = a*(1-e**2)/(1 + e*np.cos(thetas- (np.pi + psi0)))
x = np.array([r*np.cos(thetas), r*np.sin(thetas)])
plt.plot(x[0], x[1], '-.r')

def gravity(m1, m2, x):
    return -G*m1*m2/nt.norm(x)**3*x

acceleration = np.array()

for x0, y0, vx0, vy0, m in zip(x0, y0, vx0, vy0, m):
    acceleration = lambda x, t: gravity(m_star, m, x)/m
    plt.scatter(x0, y0)
    time = np.linspace(0, 0.1,  10000)
    initial_position = np.array([x0, y0])
    initial_velocity = np.array([vx0, vy0])
    xx, vv = nt.leapfrog(initial_position, initial_velocity, time, acceleration)
    plt.plot(xx[:,0], xx[:,1], linewidth = 0.8)
    plt.axis('equal')
    plt.xlabel('x'); plt.ylabel('y')
    plt.title('Planetary Orbits')
plt.show()

'''
angle = np.linspace(theta, theta + 2*np.pi, 1000)
r = a*(1-e**2)/(1 + e*np.cos(angle - (np.pi + psi)))
x = np.array([r*np.cos(angle), r*np.sin(angle)])
plt.plot(x[0], x[1], '.g')
plt.scatter(x0, y0)
'''

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

def gravity(m1, m2, x):
    #    print(G**m1*m2/nt.norm(x)**3*x)
    return -G*m1*m2/nt.norm(x)**3*x

def rotate(vector, angle):
    # Rotation matrix
    x1 = np.cos(angle)
    x2 = np.sin(angle)
    rotmatrix = np.array([[x1, -x2], [x2, x1]])
    return np.dot(vector, rotmatrix)

acc = lambda x, t: gravity(m_star, m[0], x)/m[0]

for theta, psi, a, e, x0, y0, vx0, vy0 in zip(theta0, psi0, a, e, x0, y0, vx0, vy0):
    angle = np.linspace(theta, theta + 2*np.pi, 1000)
    r = a*(1-e**2)/(1 + e*np.cos(angle))
    x = np.array([r*np.cos(angle), r*np.sin(angle)]).transpose()

    x = rotate(x, psi + np.pi)
    plt.plot(x[:, 0], x[:,1])
    plt.scatter(x0, y0)
    time = np.linspace(0, 10,  100000)
    initial_position = np.array([x0, y0])
    initial_velocity = np.array([vx0, vy0])
    xx, vv = nt.leapfrog(initial_position, initial_velocity, time, acc)
    plt.plot(xx[:,0], xx[:,1], linewidth = 0.8)

    plt.axis('equal')
    plt.xlabel('x'); plt.ylabel('y')
    plt.title('Planetary Orbits')
plt.show()

'''
This is a test for the calculated radius (from ellipse formula) vs. sqrt(x0**2 + y0**2)
Seems that it's not gonna get much closer...
Decided that the easiest way to implement an orbit of any kind was to simply
rotate an ellipse from the origin by the inclination psi.

r_calc               r_assignment
0.34088881931000203 0.3475322418676838
0.46245543984259346 0.46628905945120847
3.136171444566258 3.2266464987667134
1.9934623737098005 1.9876202428208163
0.21420545424335802 0.21626184615601463
1.1235659125452286 1.125812221107684
0.8821933895378447 0.8816067550484565
'''

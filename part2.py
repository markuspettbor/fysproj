from AST2000SolarSystem import AST2000SolarSystem
import matplotlib.pyplot as plt
import numpy as np
import numtools as nt
import prettyplot as pretty
import variables as vars

user = 'markusbp'

seed = AST2000SolarSystem.get_seed(user)
solar_system = AST2000SolarSystem(seed)


n = vars.n
x0 = vars.x0
y0 = vars.y0
vx0 = vars.vx0
vy0 = vars.vy0
a = vars.a
e = vars.e
theta0 = vars.theta0
radius = vars.radius
m_star = vars.m_star
m = vars.m
G = vars.G

def gravity(m1, m2, x):
    #    print(G**m1*m2/nt.norm(x)**3*x)
    return -G*m1*m2/nt.norm(x)**3*x

acc = lambda x, t: gravity(m_star, m[0], x)/m[0]

for a, e, theta0, x0, y0, rad, vx0, vy0 in zip(a, e, theta0, x0, y0, radius, vx0, vy0 ):
    theta = np.linspace(theta0, theta0 + 2*np.pi, 100000)
    area = np.log(2*np.pi*(rad)**2)**2/30
    r = a*(1 - e**2)/(1 + e*np.cos(theta))
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    plt.plot(x, y, linewidth = 0.8)
    plt.scatter(x0, y0, s = area)
    plt.axis('equal')

    time = np.linspace(0, 10,  20000)

    initial_position = np.array([x0, y0])
    initial_velocity = np.array([vx0, vy0])
    xx, vv = nt.leapfrog(initial_position, initial_velocity, time, acc)
    plt.plot(xx[:,0], xx[:,1], linewidth = 0.8)
    plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')

plt.title('Planetary Orbits')
plt.show()

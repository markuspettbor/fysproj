from AST2000SolarSystem import AST2000SolarSystem
import matplotlib.pyplot as plt
import numpy as np
import numtools as nt
import variables as vars

def gravity(m1, m2, x):
    return -G*m1*m2/nt.norm(x)**3*x

def kepler3(m1, m2, a):
    return np.sqrt(4*np.pi**2/(G*(m1+m2))*a**3)

def orbit(x0, v0, acc, t0, t1, steps):
    '''
    Assumes x0, v0, are arrays containing [x0, y0], [vx0, vy0]
    acc is a function defining the acceleration affecting
    the body to be simulated.
    returns x (orbital position) and v (orbital velocity)
    '''
    time = np.linspace(t0, t1, steps)
    x, v = nt.leapfrog(x0, v0, time, acc)
    return x, v, time



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

# Analytical solution
k = 100
thetas = np.array([np.linspace(theta, theta + 2*np.pi, k) for theta in theta0])
thetas = np.transpose(thetas)
r = a*(1-e**2)/(1 + e*np.cos(thetas- (np.pi + psi0)))
x_analytical = np.array([r*np.cos(thetas), r*np.sin(thetas)])

# Numerical solution

orbits = 20
stepsperorbit = 10000
period = kepler3(m_star, m, a)
t0 = 0
t1 = orbits*period
step = orbits*stepsperorbit

orbital_x = []
orbital_y = []
velocity_vx = []
velocity_vy = []
all_time = []

for m, p, index, t1 in zip(m, period, range(n), t1):
    acc = lambda r, t: gravity(m_star, m, r)/m
    period = kepler3(m_star, m, a)
    initial_x = np.array([x0[index], y0[index]])
    initial_v = np.array([vx0[index], vy0[index]])
    x, v, t = orbit(initial_x, initial_v, acc, t0, t1, step)
    orbital_x.append(x[0])
    orbital_y.append(x[1])
    velocity_vx.append(v[0])
    velocity_vy.append(v[1])
    all_time.append(t)

for x,y in zip(orbital_x, orbital_y):
    plt.plot(x, y, linewidth = 0.6)
    plt.scatter(x0, y0)
    plt.axis('equal')
plt.plot(x_analytical[0], x_analytical[1], '-.r',linewidth = 0.8)
plt.xlabel('x'); plt.ylabel('y')
plt.title('Planetary Orbits')
plt.show()

'''
acceleration = lambda x, t: gravity(m_star, m, x)/m
time = np.array([np.linspace(0, period, k*period) for period in range(1, n+1)])
initial_position = np.array([x0, y0])
initial_velocity = np.array([vx0, vy0])
xx, vv = nt.leapfrog(initial_position, initial_velocity, time, acceleration)
plt.plot(xx[:, 0], xx[:,1])
plt.show()

angle = np.linspace(theta, theta + 2*np.pi, 1000)
r = a*(1-e**2)/(1 + e*np.cos(angle - (np.pi + psi)))
x = np.array([r*np.cos(angle), r*np.sin(angle)])
plt.plot(x[0], x[1], '.g')
plt.scatter(x0, y0)

'''

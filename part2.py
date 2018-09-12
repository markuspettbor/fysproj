from AST2000SolarSystem import AST2000SolarSystem
import matplotlib.pyplot as plt
import numpy as np
import numtools as nt
import variables as vars

def gravity(m1, m2, x):
    return -G*m1*m2/nt.norm(x)**3*x

def kepler3(m1, m2, a):
    return np.sqrt(4*np.pi**2/(G*(m1+m2))*a**3)

def orbit(x0, v0, acc, t0, t1, steps, simple):
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

'''
orbits = 20
stepsperorbit = 10000
period = kepler3(m_star, m, a)
t0 = 0
t1 = orbits*period
#<<<<<<< HEAD
steps = orbits*stepsperorbit

=======
step = orbits*stepsperorbit
>>>>>>> a131a0713c5f70e7003179aedeaff5bb1763c866

orbital_x = []
orbital_y = []
velocity_vx = []
velocity_vy = []
all_time = []

<<<<<<< HEAD
'''
# First part of orbital simulation, sun at the centre of mass
'''
for m, p, index, t1, step in zip(m, period, range(n), t1, steps):
=======
for m, p, index, t1 in zip(m, period, range(n), t1):
>>>>>>> a131a0713c5f70e7003179aedeaff5bb1763c866
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
# Second part of orbital simulation, in a centre of mass reference frames

t0 = 0
M = m[3] + m_star

sun_x0 = np.array([0, 0])
sun_v0 = np.array([0, 0])
planet_x0 = np.array([x0[3], y0[3]])
planet_v0 = np.array([vx0[3], vy0[3]])

cm = m[3]/M*planet_x0 + m_star/M*sun_x0  # Centre of mass
sun_x0 = np.array([0, 0]) - cm
planet_x0 = np.array([x0[3], y0[3]]) - cm

t0 = 0
period = kepler3(m_star, m[3], a[3])
orbits = 2
stepsperorbit = 100000
t1 = orbits*period
steps = orbits*stepsperorbit
time2 = np.linspace(t0, t1, steps)
dt = time2[1] - time2[0]
mass = np.array([m_star, m[3]])

num_bodies = 2

def system_gravity(m1, m2, x):
    return -G*m1*m2/nt.norm(x, ax = 1)**3*x


def system_acceleration(m, r, index):

    mask_index = np.arange(num_bodies) != index

    r_planet = r[index]
    r_between = r[index] - r[mask_index]
    acc = system_gravity(m[mask_index], m[index], r_between)/m[index]
    #print(acc*m[index])
    acc = np.sum(acc, axis = 0)
    #print(acc)
    return acc


xx = np.zeros((len(time2), num_bodies, 2))
vv = np.zeros((len(time2), num_bodies, 2))

x = np.array([sun_x0, planet_x0])
v = np.array([sun_v0, planet_v0])


xx[0] = np.copy(x)
vv[0] = np.copy(v)

for k in range(len(time2)-1):

    for i in range(num_bodies):
        x = np.copy(xx[k])
        acc = lambda r: system_acceleration(mass, r, i)
        x[i] = xx[k,i] + vv[k,i]*dt + 0.5*acc(xx[k])*dt**2
        v[i] = vv[k,i] + 0.5*(acc(xx[k])+ acc(x))*dt

    vcm = (m[3]/M*x[1] + m_star/M*x[0]-cm)/dt
    cm = m[3]/M*x[1] + m_star/M*x[0]  # Centre of mass

    xx[k+1] = x - cm
    vv[k+1] = v

for i in range(2):
    plt.plot(xx[:,i,:][:,0],xx[:,i,:][:,1])

plt.scatter(cm[0], cm[1])
plt.axis('equal')
plt.show()



#print(m, np.sqrt(x0**2 + y0**2))

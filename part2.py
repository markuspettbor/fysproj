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

'''
# Analytical solution
k = 100
thetas = np.array([np.linspace(theta, theta + 2*np.pi, k) for theta in theta0])
thetas = np.transpose(thetas)
r = a*(1-e**2)/(1 + e*np.cos(thetas- (np.pi + psi0)))
x_analytical = np.array([r*np.cos(thetas), r*np.sin(thetas)])
plt.figure()
plt.plot(x_analytical[0], x_analytical[1])
plt.axis('equal')
#plt.show()
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

# First part of orbital simulation, sun at the centre of mass

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
# Second part of orbital simulation, in a centre of mass reference frames

# Third part of orbital simulation, in a centre of mass reference frames with more planets
def system_gravity(m1, m2, x):
    x = x.transpose()
    return (-G*m1*m2/nt.norm(x)**3*x).transpose()

def system_acceleration(m, r, index, n):
    mask_index = np.arange(n) != index
    r_planet = r[index]
    r_between = r[index] - r[mask_index]
    acc = system_gravity(m[mask_index], m[index], r_between)/m[index]

    acc = np.sum(acc, axis = 0)
    return acc

def center_of_mass(m, r):
    total_mass = np.sum(m)
    cm = np.sum(m*r.transpose(), axis = 1)/total_mass
    return cm


def n_body_problem(xx, vv, cm, vcm, mass, time, n):
    v = np.zeros(vv[0].shape)
    for k in range(len(time)-1):
        x = np.copy(xx[k])
        for i in range(n):
            acc = lambda r: system_acceleration(mass, r, i, n)
            x[i] = xx[k,i] + vv[k,i]*dt + 0.5*acc(xx[k])*dt**2
        for i in range(n):
            acc = lambda r: system_acceleration(mass, r, i, n)
            v[i] = vv[k,i] + 0.5*(acc(xx[k])+ acc(x))*dt
        cm[k+1] = center_of_mass(mass, x)
        vcm[k] = (cm[k+1]-cm[k])/dt
        xx[k+1] = x
        vv[k+1] = v
    vcm[-1] = vcm[-2] # Approximate final value
    return xx, vv, cm, vcm
'''
# One planet pulling on sun
# Init planets and suns motions
t0 = 0
mask = [3] # Selected planets
body_x0 = np.array([[0],[0]]) # Sun x0
body_v0 = np.array([[0],[0]]) # Sun v0
bodies_x0 = np.concatenate((body_x0, np.array([x0[mask], y0[mask]])), axis=1).transpose()
bodies_v0 = np.concatenate((body_v0, np.array([vx0[mask], vy0[mask]])), axis=1).transpose()

mass = np.append(m_star, m[mask])

# Time, period etc
t0 = 0
period = kepler3(mass[0], m[mask], a[mask])
orbits = 1
t1 = orbits*period    #orbits*period
stepsperorbit = 100000
steps = orbits*stepsperorbit
time2 = np.linspace(t0, t1, steps)
dt = time2[1] - time2[0]

# Shift reference frame to CM-system
cm = np.zeros((steps, 2))
cm[0] = center_of_mass(mass, bodies_x0)
vcm = np.zeros((steps, 2))
vcm[0] = np.array([0,0])
bodies_x0 = bodies_x0 - cm[0]
bodies_v0 = bodies_v0 - vcm[0]

num_bodies = len(mass)

xx = np.zeros((len(time2), num_bodies, 2))
vv = np.zeros((len(time2), num_bodies, 2))

xx[0] = np.copy(bodies_x0)
vv[0] = np.copy(bodies_v0)

xx, vv, cm = n_body_problem(xx, vv, mass, time2, num_bodies)
xx = xx.transpose()
plt.plot(xx[0,0], xx[1,0], xx[0,1], xx[1,1])
plt.scatter(cm[1:,0], cm[1:,1])
plt.axis('equal')
plt.show()

'''
# Third part of orbital simulation, in a centre of mass reference frames with more planets
mask = [0,1,2,3,4,5,6] # Selected planets
mass = np.append(m_star, m[mask])
t0 = 0
period = kepler3(mass[0], m[mask], a[mask])
stepsperorbit = 1000000
orbits = 0.01
t1 = orbits#*period*1.5
steps = int(orbits*stepsperorbit)
time2 = np.linspace(t0, t1, steps)
dt = time2[1] - time2[0]


body_x0 = np.array([[0],[0]]) # Sun x0.transpose()
body_v0 = np.array([[0],[0]]) # Sun v0
bodies_x0 = np.concatenate((body_x0, np.array([x0[mask], y0[mask]])), axis=1).transpose()
bodies_v0 = np.concatenate((body_v0, np.array([vx0[mask], vy0[mask]])), axis=1).transpose()

cm  = np.zeros((steps, 2))
vcm = np.zeros((steps, 2))
cm[0] = center_of_mass(mass, bodies_x0)
bodies_x0 = bodies_x0
bodies_v0 = bodies_v0

num_bodies = len(mass)

xx = np.zeros((len(time2), num_bodies, 2))
vv = np.zeros((len(time2), num_bodies, 2))

xx[0] = np.copy(bodies_x0)
vv[0] = np.copy(bodies_v0)
print(xx.shape)
xx, vv, cm, vcm = n_body_problem(xx, vv, cm, vcm, mass, time2, num_bodies)
xx = xx.transpose()
vv = vv.transpose()
xx_launch = np.copy(xx)
vv_launch = np.copy(vv)
print(xx.shape)
for i in range(2):
    xx[i] = xx[i] - cm[:,i]
    vv[i] = vv[i] - vcm[:,i]

xfin = xx
vfin = vv

for i in range(num_bodies):
    plt.plot(xx[0,i], xx[1,i])

plt.axis('equal')
plt.show()

def save_2Ddata(file, data):
    save = np.zeros([len(data[0])*2, len(data[0,0])])
    for i in range(len(data[0])):
        n = 2*i
        save[n] = data[0,i]
        save[n+1] = data[1,i]
    np.save(file, save.transpose())
def save_1Ddata(data, file):
    save = np.array([data, ]).transpose()
    np.save(file, save)

np.save('saved_orbits/launch_resolution/pos.npy', xx_launch)
np.save('saved_orbits/launch_resolution/vel.npy', vv_launch)
np.save('saved_orbits/launch_resolution/time.npy', time2)

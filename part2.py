from AST2000SolarSystem import AST2000SolarSystem
import matplotlib.pyplot as plt
import numpy as np
import numtools as nt

user = 'markusbp'

seed = AST2000SolarSystem.get_seed(user)
solar_system = AST2000SolarSystem(seed)

theta = np.linspace(0, 2*np.pi, 10000)

x0 = solar_system.x0
y0 = solar_system.y0
vx0 = solar_system.vx0
vy0 = solar_system.vy0

a = solar_system.a
e = solar_system.e
r0 = np.sqrt(x0**2 + y0**2)
theta0 = solar_system.omega
radius = solar_system.radius
'''
fig = plt.figure()
ax = fig.add_subplot(111, projection = 'polar')
for a, e, r0, theta0, r, in zip(a, e, r0, theta0, radius):
    area = np.log(np.pi*r**2)
    r = a*(1 - e**2)/(1 + e*np.cos(theta))
    ax.plot(theta, r)
    #ax.scatter(theta0, r0, s = area)
plt.show()
'''

m_star = solar_system.star_mass
m = solar_system.mass

m_sol = 1.989e30 #kg/solmass
mperau = 150e9 #m/AU

G = 39.478 #Nm^2/kg^2  *AUperm**2*sunmassesperthing

def gravity(m1, m2, x):
    #    print(G**m1*m2/nt.norm(x)**3*x)
    return -G*m1*m2/nt.norm(x)**3*x

t = np.linspace(0, 0.6,  10000)
acc = lambda x, t: gravity(m_star, m[0], x)/m[0]
x00 = np.array([x0[0], y0[0]])
v00 = np.array([vx0[0], vy0[0]])
x, y = nt.leapfrog(x00, v00, t, acc)

plt.plot(x[:,0], x[:,1])
plt.show()

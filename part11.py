import numpy as np
import orbit_tools as ot
import matplotlib.pyplot as plt
import variables as vars
import numtools as nt


def find_orbital_params(x_sat, x_planet, time, v_sat = 0):
    r = nt.norm(x_planet - x_sat, ax = 1)
    api  = np.max(r)
    peri = np.min(r)
    a = (api + peri)/2
    print('Periapsis:', peri)
    print('Apoapsis:', api)
    print('Semimajor:', a)
    dt = time[1] - time[0]
    rad_unit_vec = (x_planet - x_sat)/nt.norm(x_planet - x_sat)
    #print(tan_unit_vec)
    b = np.sqrt(np.min(r)*np.max(r))
    P = time[-1] - time[0]
    e = np.sqrt(1 - (b/a)**2)
    print('Period:', P)
    print('Period, Kepler 3:', ot.kepler3(vars.m[1], 0*60000/vars.solar_mass, (api + peri)/2))
    print('Eccentricity:',e)
    print('Semiminor:',b)
    print('dt:', dt)
    print('\n')


x1 = np.load('satposfororbparam.npy')[1:-2]
v1 = np.load('satvelfororbparam.npy')[1:-2]
p2 = np.load('planposfororbparam.npy')[1:-2]
torbit = np.load('timefororbparam.npy')[1:-2]

circ_point = 159

planet_radius = vars.radius[1]*1000/vars.AU_tall
# Full orbit, with ellipse and circ.
t = np.linspace(0, 2*np.pi, 100)
planet_radius = planet_radius
plt.plot(planet_radius*np.cos(t)*10**5, planet_radius*np.sin(t)*10**5, 'c', linewidth = 0.8)
plt.plot((x1[:,0] - p2[:,0])*10**5, (x1[:,1]-p2[:,1])*10**5, '--k', linewidth = 0.8)
plt.scatter((x1[:,0] - p2[:,0])[0]*10**5, (x1[:,1]-p2[:,1])[0]*10**5, c = 'k')
plt.legend(['Planet 1 Surface', 'Satellite Trajectory'], frameon = False)
plt.xlabel('x [$10^{-5}$AU]', size = 12); plt.ylabel('y [$10^{-5}$AU]', size = 12)
plt.axis('equal')
plt.show()

# Only elliptical orbit, before circ.
xell = x1[:circ_point]
pell = p2[:circ_point]
planet_radius = vars.radius[1]*1000/vars.AU_tall
find_orbital_params(xell, pell, torbit[:circ_point], v1[:circ_point])

plt.plot(planet_radius*np.cos(t)*10**5, planet_radius*np.sin(t)*10**5, linewidth = 0.8)
plt.plot((xell[:,0] - pell[:,0])*10**5, (xell[:,1]-pell[:,1])*10**5, '--k', linewidth = 0.8)
plt.legend(['Planet 1 Surface', 'Satellite Trajectory'], frameon = False)
plt.xlabel('x [$10^{-5}$AU]', size = 12); plt.ylabel('y [$10^{-5}$AU]', size = 12)
plt.axis('equal')
plt.show()

# Only final circ. orbit
xcirc = x1[circ_point + 777:]
pcirc = p2[circ_point + 777:]
find_orbital_params(xcirc, pcirc, torbit[circ_point + 777:])

plt.plot(planet_radius*np.cos(t)*10**5, planet_radius*np.sin(t)*10**5, linewidth = 0.8)
plt.plot((xcirc[:,0] - pcirc[:,0])*10**5, (xcirc[:,1]-pcirc[:,1])*10**5, '--k', linewidth = 0.8)
plt.legend(['Planet 1 Surface', 'Satellite Trajectory'], frameon = False)
plt.xlabel('x [$10^{-5}$AU]', size = 12); plt.ylabel('y [$10^{-5}$AU]', size = 12)
plt.axis('equal')
plt.show()
'''

p2= np.load('planorbforresult.npy')
x1 = np.load('optimizedorbforresult.npy')
v1 = np.load('optimizedvelforresult.npy')
opto = np.load('saved/saved_params/xopt1.npy')
optv = np.load('saved/saved_params/vopt1.npy')

x0 = np.load('unoptorbforresults.npy')
v0 = np.load('unoptvelforresults.npy')

plt.plot(x1[:,0], x1[:,1], 'k', linewidth = 0.8)
plt.plot(x0[:,0], x0[:, 1],'--', linewidth = 0.8)
plt.plot(p2[:,0], p2[:, 1], linewidth = 0.8)
plt.plot(opto.transpose()[:, 0], opto.transpose()[:,1], '-.r', linewidth = 0.8)
plt.legend(['Optimized Orbit', 'Unoptimized Orbit', 'Planet 1 Orbit', 'Optimal Orbit'], frameon = False)
plt.xlabel('x [AU]', size = 12); plt.ylabel('y [AU]', size = 12)

plt.axis('equal')
plt.show()
#t_inject = 0.5843283088441256

inject_point = np.argmin(nt.norm(p2-x1, ax = 1))
plt.plot(x1[inject_point, 0], x1[inject_point, 1],  'ok')
plt.plot(p2[inject_point, 0], p2[inject_point, 1], 'or')
plt.plot(x1[:,0], x1[:,1], 'k', linewidth = 0.8)
plt.plot(p2[:,0], p2[:, 1], linewidth = 0.8)

plt.plot(opto.transpose()[:, 0], opto.transpose()[:,1], '-.r', linewidth = 0.8)
plt.xlabel('x [AU]', size = 12); plt.ylabel('y [AU]', size = 12)

plt.legend(['Planet 1', 'Satellite' ], frameon = 0)
plt.axis('equal')
plt.show()
'''

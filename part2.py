from ast2000solarsystem import AST2000SolarSystem
import matplotlib.pyplot as plt
import numpy as np
import numtools as nt
import variables as vars
import orbit_tools as ot
from numba import jit

# First part of orbital simulation, sun at the centre of mass

def verification():
    m_star = vars.m_star
    m = vars.m
    solar_system = vars.solar_system
    theta0 = vars.theta0
    psi0 = vars.psi0
    a = vars.a
    e = vars.e
    x0 = vars.x0; y0 = vars.y0
    vx0 = vars.vx0; vy0 = vars.vy0
    # Analytical solution
    k = 10000
    thetas = np.array([np.linspace(theta, theta + 2*np.pi, k) for theta in theta0])
    thetas = np.transpose(thetas)
    r = a*(1-e**2)/(1 + e*np.cos(thetas- (np.pi + psi0)))
    x_analytical = np.array([r*np.cos(thetas), r*np.sin(thetas)])

    # Numerical solution, two body system
    orbits = 21
    stepsperorbit = 10000
    period = ot.kepler3(m_star, m, a)
    t0 = 0
    t1 = orbits*period[0]
    step = orbits*stepsperorbit
    t = np.linspace(t0, t1, step)
    x0 = np.array([[x0, y0] for x0, y0 in zip(vars.x0, vars.y0)])
    v0 = np.array([[vx0, vy0] for vx0, vy0 in zip(vars.vx0, vars.vy0)])
    x0 = np.concatenate((np.zeros((1,2)), x0))  # Set sun initial conditions
    v0 = np.concatenate((np.zeros((1,2)), v0))
    mass = np.append(m_star, vars.m)
    xx, vv = ot.patched_conic_orbits(t, mass, x0, v0)
    # Verification from MCast
    testpos = xx.transpose()[:, 1:, :]
    test = vars.solar_system
    Nyr = step/(t1)
    Tsim = t1 - t0
    time = t

    test.check_planet_positions(testpos, Tsim, Nyr)
    test.orbit_xml(testpos, time)
    plt.plot(xx[:,:,0], xx[:, :, 1], 'r')
    plt.plot(x_analytical[0], x_analytical[1], '-.r',linewidth = 0.8)
    plt.axis('equal')
    plt.xlabel('x [AU]'); plt.ylabel('y [AU]')
    plt.title('Planetary Orbits')
    plt.show()

# Third part of orbital simulation, in a centre of mass reference frames with more planets

def find_orbits():
    x0 = vars.x0
    y0 = vars.y0
    vx0 = vars.vx0
    vy0 = vars.vy0
    m_star = vars.m_star
    m = vars.m
    mask = np.arange(len(m)) # Selected planets
    mass = np.append(m_star, m[mask])
    steps = 30000
    time = np.linspace(0, 1, steps)
    body_x0 = np.array([[0],[0]]) # Sun x0
    body_v0 = np.array([[0],[0]]) # Sun v0
    _x0 = np.concatenate((body_x0, np.array([x0[mask], y0[mask]])), axis=1)
    _v0 = np.concatenate((body_v0, np.array([vx0[mask], vy0[mask]])), axis=1)
    _x0 = _x0.transpose(); _v0 = _v0.transpose()
    xx, vv, cm, vcm = ot.n_body_setup(mass, time, steps, _x0, _v0, ref_frame = 'cm')
    for i in range(len(mass)):
        plt.plot(xx[0,i], xx[1,i])
    plt.axis('equal')
    plt.show()
    # print('LENGTH TIME', len(time))
    # np.save('saved/saved_orbits/launch_resolution/pos_1cmnpy', xx)
    # np.save('saved/saved_orbits/launch_resolution/vel_1cm.npy', vv)
    # np.save('saved/saved_orbits/launch_resolution/time_1cm.npy', time)

'''
# Signal analysis
def distance_squared(x, model):
    return (x - model)**2

@jit(nopython = True)
def best_fit(t, vpec, pv, pp, pt0, best_sum, best_vr, best_p, best_t0):
    for vr in pv:
        for p in pp:
            for t0 in pt0:
                summ = np.sum((x - (vr*np.cos(2*np.pi/p*(t - t0)) + vpec))**2)
                if summ < best_sum:
                    best_sum = summ
                    best_vr = vr
                    best_p = p
                    best_t0 = t0
    return best_vr, best_p, best_t0

def smooth(signal, points):
    fit = np.ones(points)/points
    return np.convolve(signal, fit, mode = 'same')

#Exact vals:
v_pec = 0.420
vv = 0.0047396
v_r = vv
t0  = 2.5861613
P  =  3.6100834
vpp = v_r*vars.m_star/vars.m[3]

print('Exact values:')
print('Inclination: 3/7*pi')
print('Peculiar velocity: ', v_pec, '[AU/yr]')
print('Vr_star = ', v_r, '[AU/yr]')
print('P: ', P, ' [yr]')
print('t0: ', t0, '[yr]')
print('Mass: ', vars.m[3], ' [Solar masses]')
print('Planet radial velocity:  %.7f' %vpp, '[AU/yr]')
print('Radius: ', vars.radius[3], ' [km]')
print('Density: ', vars.m_normal_unit[3]/(4/3*np.pi*(vars.radius[3]*1000)**3), ' [kg/m**3]')
'''
'''
# Only for estimating number of planets
filee = np.loadtxt('saved/saved_params/rvmultiplanet.txt')
t = filee[:,0]
x = filee[:,1]
xx = smooth(x, 500)
plt.plot(t[1000:190000], xx[1000: 190000])
plt.show()
'''
'''
filee = np.loadtxt('saved/saved_params/radialvelocity.txt')
t = filee[:,0]
x = filee[:,1]


model = lambda t, vr_star, period, t0: vr_star*np.cos(2*np.pi/period*(t - t0)) + vp
steps = len(t)
vp = np.sum(x)/len(x)
v_max = (np.max(x) - np.min(x))/2
v_min = 0

p_max = np.max(t) - np.min(t)
p_min = 2*np.pi/steps

t0_max = t[-1]
t0_min = 0

stepdown = 500  # Higher stepdown means lower accuracy but greater speed

possible_t0 = np.linspace(t0_min, t0_max, steps/stepdown)
possible_p = np.linspace(p_min, p_max, steps/stepdown)
possible_v = np.linspace(v_min, v_max, steps/stepdown)

best_vr, best_p, best_t0 = best_fit(t, vp, possible_v, possible_p, possible_t0, 1e20, v_max, p_max, t0_max)
best_fit = model(t, best_vr, best_p, best_t0)
plt.plot(t, x, t, best_fit, linewidth = 0.8)
plt.xlabel('Time [h]')
plt.ylabel('Amplitude')
plt.legend(['Data', 'Best Fit for Model'], frameon = False)
plt.show()
'''

if __name__ == '__main__':
    #find_orbits()
    verification()

'''
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
'''

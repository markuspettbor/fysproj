from AST2000SolarSystem import AST2000SolarSystem
import matplotlib.pyplot as plt
import numpy as np
import numtools as nt
import variables as vars
import orbit_tools as ot

# First part of orbital simulation, sun at the centre of mass

def verification():
    m_star = vars.m_star
    m = vars.m
    solar_system = vars.solar_system
    theta0 = vars.theta0
    psi0 = vars.psi0
    a = vars.a
    e = vars.e
    # Analytical solution
    k = 1000
    thetas = np.array([np.linspace(theta, theta + 2*np.pi, k) for theta in theta0])
    thetas = np.transpose(thetas)
    r = a*(1-e**2)/(1 + e*np.cos(thetas- (np.pi + psi0)))
    x_analytical = np.array([r*np.cos(thetas), r*np.sin(thetas)])
    # Numerical solution
    orbits = 25
    stepsperorbit = 10000
    period = ot.kepler3(m_star, m, a)
    t0 = 0
    t1 = orbits*period[0]
    step = orbits*stepsperorbit
    t = np.linspace(t0, t1, step)
    orbital_x = []
    orbital_y = []
    velocity_vx = []
    velocity_vy = []
    all_time = []

    for m, index in zip(m, range(n)):
        acc = lambda r, t: gravity(m_star, m, r)/m
        period = kepler3(m_star, m, a)
        initial_x = np.array([x0[index], y0[index]])
        initial_v = np.array([vx0[index], vy0[index]])
        x, v = orbit(initial_x, initial_v, acc, t)
        orbital_x.append(x[0])
        orbital_y.append(x[1])
        velocity_vx.append(v[0])
        velocity_vy.append(v[1])
        all_time.append(t)

    testpos = np.array([orbital_x, orbital_y])
    test = vars.solar_system
    Nyr = step/(t1)
    Tsim = t1 - t0
    time = t
    test.check_planet_positions(testpos, Tsim, Nyr)
    test.orbit_xml(testpos, time)


    for x,y in zip(orbital_x, orbital_y):
        plt.plot(x, y, linewidth = 0.6)
        plt.scatter(x0, y0)
        plt.axis('equal')
    plt.plot(x_analytical[0], x_analytical[1], '-.r',linewidth = 0.8)
    plt.xlabel('x'); plt.ylabel('y')
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
    steps = 1000
    time = np.linspace(0, 1, steps)
    body_x0 = np.array([[0],[0]]) # Sun x0
    body_v0 = np.array([[0],[0]]) # Sun v0
    _x0 = np.concatenate((body_x0, np.array([x0[mask], y0[mask]])), axis=1)
    _v0 = np.concatenate((body_v0, np.array([vx0[mask], vy0[mask]])), axis=1)
    xx, vv, cm, vcm = ot.n_body_setup(mass, time, steps, _x0, _v0, ref_frame = 'cm')
    for i in range(len(mass)):
        plt.plot(xx[0,i], xx[1,i])
    plt.axis('equal')
    plt.show()

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



stepdown = 500

filee = np.loadtxt('rvmultiplanet.txt')
t = filee[:,0]
x = filee[:,1]
xx = smooth(x, 500)
plt.plot(t[1000:190000], xx[1000: 190000])
plt.show()

model = lambda t, vr_star, period, t0: vr_star*np.cos(2*np.pi/period*(t - t0)) + vp
steps = len(t)
vp = np.sum(x)/len(x)
v_max = (np.max(x) - np.min(x))/2
v_min = 0

p_max = np.max(t) - np.min(t)
p_min = 2*np.pi/steps

t0_max = t[-1]
t0_min = 0

possible_t0 = np.linspace(t0_min, t0_max, steps/stepdown)
possible_p = np.linspace(p_min, p_max, steps/stepdown)
possible_v = np.linspace(v_min, v_max, steps/stepdown)

best_vr = v_max
best_p  = p_max
best_t0 = t0_max

best_sum = np.sum(distance_squared(x, model(t, best_vr, best_p, best_t0)))

summ = best_sum

best_vr, best_p, best_t0 = best_fit(t, vp, possible_v, possible_p, possible_t0,
                                    best_sum, best_vr, best_p, best_t0)

print('Estimates:')
print('Peculiar velocity: %f [AU/yr]' %vp)
print('Vr_star: %f [AU/yr]' %best_vr)
print('P: %f [yr]' %best_p)
print('t0: %f [yr]' %best_t0)


For radialvelocity.txt:
vp = 0.28124635087721306  AU/yr
vr = 0.004498583357575979 AU/yr
P  = 100.04129857023162   yrs
t0 = 58.99871454141865    yrs


plt.plot(t, x, t, model(t, best_vr, best_p, best_t0))
plt.show()

# Lightcurve analysis section

m_s = 2.5424916122762582 # Mass of extrasolar sun

filee = np.loadtxt('rflux.txt')
tt = filee[:, 0]
xx = filee[:, 1]

plt.plot(tt, xx, tt[100:900],  smooth(xx, 10)[100:900])
plt.show()

# Simple reading of approximate plot
ttopleft = 0.05970/8766
tbotleft = 0.06682/8766
tbotright = 0.14451/8766
ttopright = 0.15120/8766

m_planet_min = m_s**(2/3)*best_vr*best_p**(1/3)/((2*np.pi*vars.G)**(1/3))
print('Minimum mass: %f [solar masses]' %m_planet_min )
v_planet = best_vr*m_s/m_planet_min
print('Planet radial velocity: %f [AU/yr] ' %v_planet)

r_planet = 0

for t1, t0 in zip((tbotleft, ttopright), (ttopleft, tbotright)):
    r_planet += (v_planet + best_vr)*(t1 - t0)/2

r_planet_avg = r_planet/2*vars.AU_tall/1000 # Convert to km
print('Radius: %f [km] or %f [km]'  %(r_planet_avg, r_planet_avg*1000))

rho_est = m_planet_min*1.989e30/(4/3*np.pi*(r_planet_avg*1000)**3) # kg/m**3
print('Density: %f [kg/m**3] or %f [kg/m**3]'  %(rho_est, rho_est/1000**3))
'''


if __name__ == '__main__':
    find_orbits()
    #verification()

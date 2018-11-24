from ast2000solarsystem import AST2000SolarSystem
import matplotlib.pyplot as plt
import numpy as np
import numtools as nt
import variables as vars
import orbit_tools as ot
from numba import jit

m_star = vars.m_star
m_sat = vars.satellite
m = vars.m
solar_system = vars.solar_system
theta0 = vars.theta0
psi0 = vars.psi0
a = vars.a
e = vars.e
x0 = vars.x0; y0 = vars.y0
vx0 = vars.vx0; vy0 = vars.vy0


def kepler3(m1, m2, a):
    return np.sqrt(4*np.pi**2/(vars.G*(m1+m2))*a**3)

def keplercheck():
    #compare two areas one where close to aphelion, one where close to perihelion
    #INPUTS: positions_vect, velocities_vect, time_vect..
    a = vars.a
    m_star = vars.m_star
    m = vars.m
    P = ot.kepler3(m_star, m, a)
    t = np.load('saved/saved_orbits/time_verif.npy')
    x = np.load('saved/saved_orbits/pos_verif.npy')
    v = np.load('saved/saved_orbits/vel_verif.npy')
    dA = np.zeros(len(x[0]))#, len(x[0][0])])
    dt = t[1]-t[0]
    #print(dt)
    for n in range(int(len(x[0])) - 1):
        vel = v[2:, n + 1, :].transpose() #9999, 8, 2
        pos = x[2:, n + 1, :].transpose()
        print('Planet %i:' %(n))
        dist = nt.norm(pos)
        apoapsis = np.argmax(dist)
        periapsis = np.argmin(dist)
        buel_api = nt.norm(pos[:,apoapsis+1] - pos[:,apoapsis-1])
        area_api = dist[apoapsis]*buel_api/2
        area_api_v2 = 1/2*dist[apoapsis]**2*nt.angle_between(pos[:,apoapsis], pos[:,apoapsis+1])
        buel_peri = nt.norm(pos[:,periapsis+1] - pos[:,periapsis-1])
        area_peri = dist[periapsis]*buel_peri/2
        area_peri_v2 = 1/2*dist[periapsis]**2*nt.angle_between(pos[:,periapsis], pos[:,periapsis+1])
        print('Apoapsis - Periapsis =', area_api - area_peri) #larger numerical erros in v2 du to small angles i think
        print('Apoapsis / Periapsis =', area_api / area_peri) #larger numerical erros in v2 du to small angles i think
        print('Distance traveled apoapsis:', buel_api)
        print('Distance traveled periapsis:', buel_peri)
        print('Mean speed apoapsis:', buel_api/(2*dt))
        print('Mean speed periapsis:', buel_peri/(2*dt))
        plt.plot(t[2:], dist)
        print('')

    plt.show()

    #3rd law
    P_k3 = np.sqrt(a**3)
    indeks = np.argmin(P_k3/P)
    comparison = P_k3/P/(P_k3[indeks]/P[indeks])
    print(comparison)
    print(max(comparison-1)/max(comparison)*100,'%')


    '''
    print(apoapsis[0])
    #measure area perihelion
    for i in num:
        rad = dist[i][apoapsis[i]]
        dx = nt.norm(x[i][apoapsis[i]+1]- x[i][apoapsis[i]])
        dA[i] = rad*dx/2
    print(dA[0])
    for i in num:
        rad = dist[i][periapsis[i]]
        dx = abs(dist[i][periapsis[i]+1]- dist[i][periapsis[i]])
        dA[i] = rad*dx/2
    print(dA[0])#measure area, dA =
    '''


def verification():
    # First part of orbital simulation, sun at origin
    # Analytical solution
    k = 10000 * 21
    thetas = np.array([np.linspace(theta, theta + 2*np.pi, k) for theta in theta0])
    thetas = np.transpose(thetas)
    r = a*(1-e**2)/(1 + e*np.cos(thetas- (np.pi + psi0)))
    x_analytical = np.array([r*np.cos(thetas), r*np.sin(thetas)])

    # Numerical solution, two body system
    orbits = 31
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
    #np.save('saved/saved_orbits/pos_verif.npy', xx)
    #np.save('saved/saved_orbits/vel_verif.npy', vv)
    #np.save('saved/saved_orbits/time_verif.npy', t)
    #keplercheck(t,xx,vv)
    #test.check_planet_positions(testpos, Tsim, Nyr)
    #test.orbit_xml(testpos, time)
    plt.plot(xx[:,:,0], xx[:, :, 1], 'k', linewidth = 0.8)
    plt.title('Planetary Orbits, Numerical Solution')
    plt.xlabel('x [AU]', size = 12); plt.ylabel('y [AU]', size = 12)
    plt.axis('equal')

    plt.figure()
    plt.plot(x_analytical[0], x_analytical[1], '--k',linewidth = 0.8)
    plt.title('Planetary Orbits, Analytical Solution')
    plt.axis('equal')
    plt.xlabel('x [AU]', size = 12); plt.ylabel('y [AU]', size = 12)
    plt.show()

def find_orbits():
    # Third part of orbital simulation, in a centre of mass reference frames with more planets
    mask = np.arange(len(m)) # Selected planets
    mass = np.append(m_star, m[mask])
    period  = ot.kepler3(m_star, m, a)[0]
    orbits = 21
    stepsperorbit = 10000
    t1 = orbits*period
    steps = orbits*stepsperorbit
    time = np.linspace(0, t1, steps)
    body_x0 = np.array([[0],[0]]) # Sun x0
    body_v0 = np.array([[0],[0]]) # Sun v0
    _x0 = np.concatenate((body_x0, np.array([x0[mask], y0[mask]])), axis=1)
    _v0 = np.concatenate((body_v0, np.array([vx0[mask], vy0[mask]])), axis=1)
    _x0 = _x0.transpose(); _v0 = _v0.transpose()
    xx, vv, cm, vcm = ot.n_body_setup(mass, time, steps, _x0, _v0, ref_frame = 'cm')
    # Find orbits for n-body system (much more fun)
    #for i in range(len(mass)):
    #    plt.plot(xx[0,i], xx[1,i])

    # Find two-body system for just star and planet 3
    mask = np.array([0, 3])
    mass2 = mass[mask]
    period  = ot.kepler3(m_star, m, a)[3]
    orbits = 21
    stepsperorbit = 10000
    t1 = orbits*period
    steps = orbits*stepsperorbit
    time = np.linspace(0, t1, steps)
    x02 = _x0[mask]
    v02 = _v0[mask]
    x2, v2, cm, vcm = ot.n_body_setup(mass2, time, steps, x02, v02, ref_frame = 'cm')
    r = nt.norm(x2[:, 1] - x2[:, 0])
    v = nt.norm(v2[:, 1] - v2[:, 0])
    total_energy = ot.energy_cm(m_star, mass2[1], v, r )
    max_dev_energy = np.max(np.abs(total_energy)) - np.min(np.abs(total_energy))
    print('Deviation, energy in cm system:', max_dev_energy)
    print('Average energy:', np.sum(total_energy)/len(total_energy))
    for i in range(2):
        plt.plot(x2[0, i], x2[1, i], '--k')
    plt.axis('equal')
    plt.show()
    # print('LENGTH TIME', len(time))
    # np.save('saved/saved_orbits/launch_resolution/pos_1cmnpy', xx)
    # np.save('saved/saved_orbits/launch_resolution/vel_1cm.npy', vv)
    # np.save('saved/saved_orbits/launch_resolution/time_1cm.npy', time)

#EXTRATERRESTRIALS
def radial_velocity_curve(vel_sun):
    inc = 3/7*np.pi
    vel_pec = 0.420
    vel_radial_data = nt.create_radial_velocity(vel_sun, vel_pec, inc)
    return vel_radial_data

def create_light_curve(pos, rp, rs): #position, time, radius planet, radius sun
    x = pos[0]
    y = pos[1]
    area_planet = np.pi*rp**2
    area = np.zeros(len(x))
    n = np.linspace(0,len(x)-1, len(x), dtype=int)
    for x, y, n in zip(x, y, n):
        if y < 0:
            area[n] = 0
        elif abs(x) >= rs+rp:
            area[n] = 0
        elif abs(x) <= rs-rp:
            area[n] = area_planet
        elif nt.norm([x,rp]) > rs:
            if x < 0: #innenfor venstre
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                + c*np.sqrt(rs**2-c**2) + rs**2*np.arcsin(c/rs) + 2*c*x
                area[n] = areal
            else:   #innenfor høyre
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                + c*np.sqrt(rs**2-c**2) + rs**2*np.arcsin(c/rs) - 2*c*x
                area[n] = areal
        else:   #utenfor venstre
            if x < 0:
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal_inv = +c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                - c*np.sqrt(rs**2-c**2) - rs**2*np.arcsin(c/rs) - 2*c*x
                area[n] = area_planet - areal_inv
            else: #utenfor høyre
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal_inv = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                - c*np.sqrt(rs**2-c**2) - rs**2*np.arcsin(c/rs) + 2*c*x
                area[n] = area_planet - areal_inv
    return area

def radial_velocity_function():
    save_orbits = False

    if save_orbits:
        print('save start')
        mask = np.arange(len(m)) # Selected planets
        mass = np.append(m_star, m[mask])
        period  = ot.kepler3(m_star, m, a)[0]
        t1 = ot.kepler3(m_star, m, a)[3]
        stepsperorbit = 10000
        steps = int(stepsperorbit*t1)
        time = np.linspace(0, t1, steps)
        body_x0 = np.array([[0],[0]]) # Sun x0
        body_v0 = np.array([[0],[0]]) # Sun v0
        _x0 = np.concatenate((body_x0, np.array([x0[mask], y0[mask]])), axis=1)
        _v0 = np.concatenate((body_v0, np.array([vx0[mask], vy0[mask]])), axis=1)
        _x0 = _x0.transpose(); _v0 = _v0.transpose()
        xx, vv, cm, vcm = ot.n_body_setup(mass, time, steps, _x0, _v0, ref_frame = 'cm')
        np.save('saved/saved_orbits/extraterrestrials/time.npy', time)
        np.save('saved/saved_orbits/extraterrestrials/xx.npy', xx)
        np.save('saved/saved_orbits/extraterrestrials/vv.npy', vv)
        np.save('saved/saved_orbits/extraterrestrials/cm.npy', cm)
        np.save('saved/saved_orbits/extraterrestrials/vcm.npy', vcm)
        print('save stop')

    time = np.load('saved/saved_orbits/extraterrestrials/time.npy')
    xx = np.load('saved/saved_orbits/extraterrestrials/xx.npy')
    vv = np.load('saved/saved_orbits/extraterrestrials/vv.npy')
    cm = np.load('saved/saved_orbits/extraterrestrials/cm.npy')
    vcm = np.load('saved/saved_orbits/extraterrestrials/vcm.npy')
    plt.rcParams.update({'font.size': 18})
    vel_radial_data = radial_velocity_curve(vv[:,0])
    plt.plot(time, vel_radial_data, 'k', linewidth = 0.1)
    plt.xlabel('Time [years]')
    plt.ylabel('Radial Velocity [AU/year]')
    plt.show()

    #def save_data():
        #velocity_radial_data_sending = np.array([time, vel_radial_data]).transpose()
        #np.savetxt('velocity_radial_one_planet_v3.txt', velocity_radial_data_sending)

        #flux_relative_data_sending = np.array([time[start:stop], flux_relative_data[start:stop]]).transpose()
        #np.savetxt('flux_relative_V2.txt', flux_relative_data_sending)
    #save_data()

def light_curve_function():
    save_orbits = False

    if save_orbits:
        print('save start')
        mask = np.arange(len(m)) # Selected planets
        mass = np.append(m_star, m[mask])
        period  = ot.kepler3(m_star, m, a)[0]
        orbits = 1
        stepsperorbit = 1000000
        t1 = orbits*period
        steps = int(orbits*stepsperorbit)
        time = np.linspace(0, t1, steps)
        body_x0 = np.array([[0],[0]]) # Sun x0
        body_v0 = np.array([[0],[0]]) # Sun v0
        _x0 = np.concatenate((body_x0, np.array([x0[mask], y0[mask]])), axis=1)
        _v0 = np.concatenate((body_v0, np.array([vx0[mask], vy0[mask]])), axis=1)
        _x0 = _x0.transpose(); _v0 = _v0.transpose()
        xx, vv, cm, vcm = ot.n_body_setup(mass, time, steps, _x0, _v0, ref_frame = 'cm')
        np.save('saved/saved_orbits/extraterrestrials/time_short.npy', time)
        np.save('saved/saved_orbits/extraterrestrials/xx_short.npy', xx)
        np.save('saved/saved_orbits/extraterrestrials/vv_short.npy', vv)
        np.save('saved/saved_orbits/extraterrestrials/cm_short.npy', cm)
        np.save('saved/saved_orbits/extraterrestrials/vcm_short.npy', vcm)
        print('save stop')

    time = np.load('saved/saved_orbits/extraterrestrials/time_short.npy')
    xx = np.load('saved/saved_orbits/extraterrestrials/xx_short.npy')
    vv = np.load('saved/saved_orbits/extraterrestrials/vv_short.npy')
    cm = np.load('saved/saved_orbits/extraterrestrials/cm_short.npy')
    vcm = np.load('saved/saved_orbits/extraterrestrials/vcm_short.npy')


    m1 = vars.m
    m2 = vars.m_star
    orbit = np.sqrt(vars.x0**2 + vars.y0**2)
    radius_planets = vars.radius*1000/vars.AU_tall
    radius_sun = vars.radius_star*1000/vars.AU_tall
    index = np.argmax(m1*m2/(orbit**2))
    radius_heavy = radius_planets[index]
    mask = [0]#,1,2,3,4,5,6]#....
    pos_planets = np.copy(xx[:,1:])
    area_sun = np.pi*radius_sun**2
    area_covered = np.zeros(len(xx[0,0]))
    area_covered_saved = []
    for mask_index in mask:
        area = create_light_curve(pos_planets[:,mask_index], radius_planets[mask_index], radius_sun)
        area_covered += area
        area_covered_saved.append(area)
    flux_relative = (area_sun-area_covered)/area_sun
    mu = 0
    sig = max(area_covered/area_sun*0.2)
    #sig = 0.01
    flux_relative_data = nt.noiceify(flux_relative, mu, sig)
    start = int(len(flux_relative_data)*0)#*(1.681/2))
    stop = int(len(flux_relative_data))#*(1.686/2))
    plt.plot(time[start:stop], flux_relative[start:stop], 'r')
    plt.plot(time[start:stop], flux_relative_data[start:stop], '-b', linewidth = 0.1)
    plt.legend(['model', 'measurement'])
    plt.xlabel('time [AU]')
    plt.ylabel('flux')
    plt.show()

    #def save_data():
        #velocity_radial_data_sending = np.array([time, vel_radial_data]).transpose()
        #np.savetxt('velocity_radial_one_planet_v3.txt', velocity_radial_data_sending)

        #flux_relative_data_sending = np.array([time[start:stop], flux_relative_data[start:stop]]).transpose()
        #np.savetxt('flux_relative_V2.txt', flux_relative_data_sending)
    #save_data()


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
    #verification()
    #keplercheck()
    #radial_velocity_function()
    #light_curve_function()
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

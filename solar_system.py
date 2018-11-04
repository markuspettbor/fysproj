import numpy as np
import numtools as nt
import orbit_tools as ot

class SolarSystem:
    def __init__(self):
        self.planets = []
        self.suns = []
        self.satellites = []
        self.bodies = []

    def addPlanet(self, planet):
        self.planets.append(planet)
        self.bodies.append(planet)

    def addSun(self, sun):
        self.suns.append(sun)
        self.bodies.append(sun)

    def addSatellite(self, satellite):
        self.satellites.append(satellite)
        self.bodies.append(satellite)

    def find_orbits(self, time, ref_frame, mass, x0, v0):
        steps = len(time)
        return ot.n_body_setup(mass, time, steps, x0, v0, ref_frame)

    def find_two_body_orbits(self, time, mass, x0, v0, ref_frame = 0):
        '''
        time is time array
        ref_frame is the index of the planet/star of the desired reference
        frame, default is sun (index 0).
        Assumes mass is array of masses, with index 0 corresponding to sun
        x0 is array of initial positions, assumed to be of shape
        (1, number of planets, 2).
        v0 is array of initial velocities, same shape.
        '''
        steps = len(time)

        x = np.zeros((steps, x0.shape[0], x0.shape[1]))
        v = np.zeros(x.shape)
        for i in range(1, len(mass)):
            acc = lambda r, t: ot.gravity(mass[0], mass[i], r)/mass[i]
            a, b = ot.orbit(x0[i], v0[i], acc, time)
            x[:, i], v[:, i] = a.transpose(), b.transpose()
            x[:, i] = x[:,i] - x[:,ref_frame]
            v[:, i] = v[:,i] - v[:,ref_frame]
        return x, v


class Planet:
    def __init__(self, mass, radius, position, velocity, name):
        self.mass = mass
        self.radius = radius
        self.position = position
        self.velocity = velocity
        self.name = name

class Sun:
    def __init__(self, mass, radius, position, velocity, temperature, name):
        self.mass = mass
        self.radius = radius
        self.position = position
        self.velocity = velocity
        self.temperature = temperature
        self.name = name

class Satellite:
    def __init__(self, mass, position, velocity, name):
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.name = name

    def addLaunch_codes(self, boost_time, delta_v):
        self.boost_time = boost_time
        self.delta_v = delta_v

if __name__ == '__main__':
    import variables as vars
    m_star = vars.m_star
    m = vars.m
    x0 = vars.x0
    y0 = vars.y0
    vx0 = vars.vx0
    vy0 = vars.vy0
    radius = vars.radius
    m_sat = vars.satellite/vars.solar_mass

    names = ['Dum','og', 'Deilig', 'Juba', 'juba', 'Pizzatryne', 'Verdens ende']
    sol = SolarSystem()
    sun = Sun(m_star, 0.001, np.array([0,0]), np.array([0,0]), 1000, 'Sol')
    sol.addSun(sun)

    for name, xx0, yy0, vvx0, vvy0, mass, r in zip(names, x0, y0, vx0, vy0, m, radius):
        x = np.array([xx0, yy0])
        v = np.array([vvx0, vvy0])
        planet = Planet(mass, r, x, v, name)
        sol.addPlanet(planet)

    steps = 30000
    t3 = 0.6
    tol = 0.0001
    time = np.linspace(0, t3, steps)
    import matplotlib.pyplot as plt

    # Generate super-orbits, hopefully only once.
    print('Simulating orbits.')
    mass = np.array([body.mass for body in sol.bodies])
    x0 = np.array([body.position for body in sol.bodies])
    v0 = np.array([body.velocity for body in sol.bodies])
    #xx, vv = sol.find_two_body_orbits(time, mass, x0, v0)
    #np.save('saved/saved_params/xx_sol1.npy', xx)
    #np.save('saved/saved_params/vv_sol1.npy', vv)
    print('Done')
    xx = np.load('saved/saved_params/xx_sol1.npy')
    vv = np.load('saved/saved_params/vv_sol1.npy')

    m_t = np.append(mass, m_sat)

    dw, tw, t_cept, semimajor = ot.trajectory(m_t, xx, vv, 1, -1, 2, 0, time, False, tol)
    print('Launch Parameters:', dw[-1], tw, t_cept[-1])

    # Simulating launch region with higher res.
    t_launch = tw[-1]
    steps2 = 10000
    launch_duration = 0.0001

    t1 = t_launch - t_launch/steps - launch_duration
    t2 = t_launch + launch_duration

    phase1 = np.linspace(0, t1, steps)
    phase2 = np.linspace(t_launch - launch_duration, t2 - (t2-t_launch)/steps2, steps2)
    phase3 = np.linspace(t2, t3, steps)
    full_time = np.concatenate((phase1, phase2, phase3))

    #xx, vv = sol.find_two_body_orbits(full_time, mass, x0, v0)
    #np.save('saved/saved_params/xx_sol2.npy', xx)
    #np.save('saved/saved_params/vv_sol2.npy', vv)

    xx = np.load('saved/saved_params/xx_sol2.npy')
    vv = np.load('saved/saved_params/vv_sol2.npy')

    dv = np.zeros((len(full_time), 2))
    t_launch_indx = len(phase1)
    cept = np.argmin(np.abs(full_time - t_cept[-1]))
    dv[t_launch_indx] = dw[-1]*nt.unit_vector(vv[t_launch_indx, 1])
    print('Go for launch')
    tol2 = 1e-6
    dw, tw, t_cept, semimajor = ot.trajectory(m_t, xx, vv, 1, -1, 2, 0, full_time, False, tol2)

    import launch
    launch.test(testing = True)
    xx = xx.transpose()
    vv = vv.transpose()

    xx = xx.transpose()
    vv = vv.transpose()
    print(dw[-1], tw, t_cept[-1])
    print('Initial launch index:', t_launch_indx)
    t_launch_indx = np.argmin(abs(full_time - tw[-1]))
    print('New launch index:', t_launch_indx)
    dv[t_launch_indx] = dw[-1]*nt.unit_vector(vv[t_launch_indx, 1])
    cept = np.argmin(abs(full_time - t_cept[-1]))
    print('Index of intercept:',cept)

    def find_launch_sequence(maxiter = 5):
        intercept = cept - t_launch_indx
        host = 1
        time2 = full_time[t_launch_indx:]
        planet_x = xx[t_launch_indx:]
        planet_v = vv[t_launch_indx:]
        dvv = dv[t_launch_indx:]
        x0_sat = planet_x[0, host] + vars.radius[0]*1000/vars.AU_tall*nt.unit_vector(planet_v[0, host])
        v0_sat = planet_v[0, host]
        xv = ot.vis_viva(mass[0], nt.norm(x0_sat), semimajor[-1])

        deviation = np.pi/2-nt.angle_between(v0_sat, x0_sat)

        dvv[0] = nt.rotate(dvv[0], deviation)

        v00_sat = xv*nt.rotate(nt.unit_vector(v0_sat), deviation)

        acc_optimal = lambda r, t: ot.gravity(mass[0], m_sat, r)/m_sat
        opt_orb, opt_vel = ot.orbit(x0_sat, v00_sat, acc_optimal, time2)
        opt_vel = opt_vel.transpose()
        opt_orb = opt_orb.transpose()
        x_sat, v_sat, dv_used = ot.n_body_sat(planet_x, mass, time2, dvv, x0_sat, v0_sat, m_sat)
        print(v_sat)
        print('Closest approach:', min(nt.norm(planet_x[:, 2] - x_sat, ax = 1)))
        opt_dev = nt.norm(x_sat - opt_orb.transpose(), ax = 1)
        boost_thresh = ot.grav_influence(mass[0], mass[1], x_sat, 1)
        dist_to_host = nt.norm(planet_x[:, 1] - x_sat, ax = 1)
        boost_time = np.where(dist_to_host > boost_thresh)[0][0]
        print('Time of boost:', time2[boost_time])
        x_sat, v_sat, dv_used = ot.n_body_sat(planet_x, mass, time2, dvv, x0_sat, v0_sat, m_sat, opt_vel, opt_orb, time2[boost_time], True)
        print('Closest approach:',min(nt.norm(planet_x[:, 2] - x_sat, ax = 1)))
        print('Optimal closest approach:', ot.grav_influence(mass[0], mass[2], x_sat)[intercept])
        closest  = np.argmin(nt.norm(planet_x[:,2]- x_sat, ax=1))
        print('total dv:', np.sum(nt.norm(dv_used)))
        dv_used[closest-3000:] = dv_used[closest-3000:]*0
        x_sat, v_sat, dv_used = ot.n_body_sat(planet_x, mass, time2, dv_used, x0_sat, v0_sat, m_sat, opt_vel, opt_orb, time2[boost_time], False)
        print('Closest approach:',min(nt.norm(planet_x[:, 2] - x_sat, ax = 1)))
        print('total dv:', np.sum(nt.norm(dv_used)))
        closest  = np.argmin(nt.norm(planet_x[:,2]- x_sat, ax=1))
        return x_sat, v_sat, dvv, opt_orb, closest, dv_used

    xss, xvv, dv, opt_orb, close, dv_used = find_launch_sequence(5)
    closest = t_launch_indx + close

    time3 = np.linspace(t3, 0.605, 30000)
    x0 = xx[closest]
    v0 = vv[closest]

    #xx, vv = sol.find_two_body_orbits(time3, mass, x0, v0)
    #np.save('saved/saved_params/xx_sol3.npy', xx)
    #np.save('saved/saved_params/vv_sol3.npy', vv)
    xx = np.load('saved/saved_params/xx_sol3.npy')
    vv = np.load('saved/saved_params/vv_sol3.npy')
    x0_sat = xss[close]
    v0_sat = xvv[close]
    dv3 = np.zeros((len(time3), 2))
    dir = nt.unit_vector(vv[0, 2])
    dir2 = nt.unit_vector(xx[0,2] - x0_sat)
    dir2 = nt.rotate(dir2, np.pi/2)
    dv3[0] =  -v0_sat + vv[0,2] - ot.vis_viva(mass[2], nt.norm(x0_sat-xx[0,2]), nt.norm(x0_sat-xx[0,2]))*dir2



    xss, v_sat, dv_used = ot.n_body_sat(xx, mass, time3, dv3, x0_sat, v0_sat, m_sat)

    plt.plot(xss[:,0] - xx[:, 2,0], xss[:,1]-xx[:, 2, 1], c = 'k')
    plt.axis('equal')
    plt.show()

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

    def update_orbits(self, dt):
        for body in self.bodies:
            acc = lambda r: ot.system_acceleration(masses, r, i, n)
            x[i] = xx[k,i] + vv[k,i]*dt + 0.5*acc(xx[k])*dt**2
        for i in range(n):
            acc = lambda r: ot.system_acceleration(masses, r, i, n)
            v[i] = vv[k,i] + 0.5*(acc(xx[k])+ acc(x))*dt

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

    steps = 10000
    t3 = 0.6
    tol = 0.00001
    time = np.linspace(0, t3, steps)
    import matplotlib.pyplot as plt

    # Generate super-orbits, hopefully only once.
    print('Simulating orbits.')
    mass = np.array([body.mass for body in sol.bodies])
    x0 = np.array([body.position for body in sol.bodies])
    v0 = np.array([body.velocity for body in sol.bodies])
    '''
    xx, vv, nada, zipp = sol.find_orbits(time, 'sol', mass, x0, v0)

    np.save('saved/saved_params/xx_sol1.npy', xx)
    np.save('saved/saved_params/vv_sol1.npy', vv)
    print('Done')
    '''
    xx = np.load('saved/saved_params/xx_sol1.npy')
    vv = np.load('saved/saved_params/vv_sol1.npy')

    xx = xx.transpose()
    vv = vv.transpose()
    m_t = np.append(mass, m_sat)

    dw, tw, t_cept, semimajor = ot.trajectory(m_t, xx, vv, 1, -1, 2, 0, time, False, tol)
    print('Launch Parameters:', dw[-1], tw[-1], t_cept[-1])

    # Simulating launch region with higher res.
    t_launch = tw[0]
    steps2 = 10000
    launch_duration = 0.001

    t1 = t_launch - t_launch/steps - launch_duration
    t2 = t_launch + launch_duration

    phase1 = np.linspace(0, t1, steps)
    phase2 = np.linspace(t_launch - launch_duration, t2 - (t2-t_launch)/steps2, steps2)
    phase3 = np.linspace(t2, t3, steps)
    full_time = np.concatenate((phase1, phase2, phase3))
    '''
    xx, vv, nada, zipp = sol.find_orbits(full_time, 'sol', mass, x0, v0)

    np.save('saved/saved_params/xx_sol2.npy', xx)
    np.save('saved/saved_params/vv_sol2.npy', vv)
    '''
    xx = np.load('saved/saved_params/xx_sol2.npy')
    vv = np.load('saved/saved_params/vv_sol2.npy')

    xx = xx.transpose()
    vv = vv.transpose()
    dv = np.zeros((len(full_time), 2))
    t_launch_indx = len(phase1)
    cept = np.argmin(np.abs(full_time - t_cept[-1]))
    dv[t_launch_indx] = dw[-1]*nt.unit_vector(vv[t_launch_indx, 1])
    print('Go for launch')
    tol2 = 1e-6
    dw, tw, t_cept, semimajor = ot.trajectory(m_t, xx, vv, 1, -1, 2, 0, full_time, False, tol2)
    print(dw[-1], tw[-1], t_cept[-1])
    print(t_launch_indx)
    t_launch_indx = np.argmin(abs(full_time - tw[-1]))
    print(t_launch_indx)
    dv[t_launch_indx] = dw[-1]*nt.unit_vector(vv[t_launch_indx, 1])

    cept = np.argmin(abs(full_time - t_cept[-1]))
    print(t_cept)
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
        print(np.arctan2(planet_x[intercept, 2, 1], planet_x[intercept, 2, 0]) - np.arctan2(opt_orb[1,0], opt_orb[0, 0]))

        print(np.arctan2(opt_orb[1,intercept], opt_orb[0, intercept])-np.arctan2(opt_orb[1,0], opt_orb[0, 0]))


        boost = 0.1
        x_sat = np.zeros((len(time2), 2))
        v_sat = np.zeros((len(time2), 2))
        x_sat[0] = x0_sat
        v_sat[0] = v0_sat
        opt_vel = opt_vel.transpose()

        def acc(r_sat, r):
            r_between = r_sat - r
            rr1 = nt.norm(r_between, ax = 1)
            acc = 0
            for mm, rr, rb in zip(mass, rr1, r_between):
                acc1 = -vars.G*mm/rr**3*rb
                acc += acc1
            return acc

        for k in range(len(time2) -1):
            dt = time2[k+1] - time2[k]
            acc1 = acc(x_sat[k], planet_x[k])
            diff = v_sat[k] - opt_vel[k]
            if k < intercept:
                v_sat[k] = v_sat[k] #- diff

            x_sat[k+1] = x_sat[k] + v_sat[k]*dt + 0.5*acc1*dt**2
            acc2 = acc(x_sat[k+1], planet_x[k+1])
            v_sat[k+1] = v_sat[k] + 0.5*(acc1 + acc2)*dt + dvv[k]
        print(nt.norm(planet_x[intercept, 2] - opt_orb[:, intercept]))
        print(nt.norm(planet_x[intercept, 2] - x_sat[intercept]))
        print(nt.norm(x_sat - opt_orb.transpose(), ax = 1) < 0.001)
        return x_sat, dvv, opt_orb
    xss, dv, opt_orb = find_launch_sequence(5)




    a = nt.norm(xx[cept,2]-xss[cept - t_launch_indx])

    xss = xss.transpose()
    for i in range(8):
        plt.plot(xx[:,i,0], xx[:, i, 1])
        plt.axis('equal')
    plt.plot(xss[0], xss[1], c = 'k')
    plt.scatter(opt_orb[0, cept - t_launch_indx], opt_orb[1, cept-t_launch_indx])
    plt.plot(opt_orb[0], opt_orb[1], c = 'm')
    plt.scatter(xss[0, cept - t_launch_indx], xss[1, cept - t_launch_indx], c = 'r')
    plt.scatter(xx[cept, 2, 0], xx[cept, 2, 1], c = 'k')
    plt.show()

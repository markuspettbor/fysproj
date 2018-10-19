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
    t1 = 0.59
    tol = 0.0005
    time = np.linspace(0, t1, steps)
    import matplotlib.pyplot as plt

    # Generate super-orbits, hopefully only once.
    print('Simulating orbits.')
    mass = np.array([body.mass for body in sol.bodies])
    x0 = np.array([body.position for body in sol.bodies])
    v0 = np.array([body.velocity for body in sol.bodies])
    xx, vv, nada, zipp = sol.find_orbits(time, 'sol', mass, x0, v0)
    '''
    np.save('saved/saved_params/xx_1mill.npy', xx)
    np.save('saved/saved_params/vv_1mill.npy', vv)
    '''
    print('Done')

    #xx = np.load('saved/saved_params/xx_1mill.npy')
    #vv = np.load('saved/saved_params/vv_1mill.npy')
    xx = xx.transpose()
    vv = vv.transpose()
    m_t = np.append(mass, m_sat)
    try:
        dw, tw, t_cept = ot.trajectory(m_t, xx, vv, steps, 1, -1, 2, 0, time, False, tol, 10000)
        print('Launch Parameters:', dw[0], tw[0], t_cept)
    except IndexError:
        print('No launch window found')

    # Simulating launch region with higher res.
    t_launch = tw[0]
    steps2 = 10000
    launc_duration = 0.001

    t1 = t_launch - t_launch/steps
    t2 = t_launch + launc_duration
    t3 = 0.60

    phase1 = np.linspace(0, t1, steps)
    phase2 = np.linspace(t_launch, t2 - t2/steps2, steps2)
    phase3 = np.linspace(t2, t3, steps)
    full_time = np.concatenate((phase1, phase2, phase3))
    xx, vv, nada, zipp = sol.find_orbits(full_time, 'sol', mass, x0, v0)
    x0_sat = x0[1] + vars.radius[0]*1000/vars.AU_tall
    v0_sat = v0[1]
    xx = xx.transpose()
    vv = vv.transpose()
    dv = np.zeros((len(full_time), 2))
    t_launch_indx = len(phase1) + 1
    cept = np.argwhere(abs(t_cept - full_time) < 0.001)[0,0]
    dv[t_launch_indx] = dw[0]*nt.unit_vector(vv[t_launch_indx, 1])
    tol_dist = 0.0005

    print('Go for launch')
    def find_launch_sequence(maxiter = 5):
        intercept = cept - t_launch_indx
        host = 1
        time = full_time[t_launch_indx:]
        planet_x = xx[t_launch_indx:]
        planet_v = vv[t_launch_indx:]
        dvv = dv[t_launch_indx:]
        x0_sat = planet_x[0, host] + vars.radius[0]*1000/vars.AU_tall*nt.unit_vector(planet_v[0, host])
        v0_sat = planet_v[0, host]
        xs, vs = ot.n_body_sat(planet_x, mass, time, dvv, x0_sat, v0_sat, m_sat)
        for t in range(5):
            dvv = dv[t_launch_indx:]
            boost = -0.01
            for i in range(maxiter):
                xs, vs = ot.n_body_sat(planet_x, mass, time, dvv, x0_sat, v0_sat, m_sat)
                d = dvv[0] + boost
                dvv[0] = nt.rotate(d, -2*np.pi/360*t)
                dist = nt.norm(planet_x[intercept, 2]- xs[intercept])
                print(dist)
        return xs, dvv
    '''
    def bazbaz(maxiter = 20):
        xs, vs = ot.n_body_custom(mass, full_time, xx, vv, 1, dv, False, x0_sat, v0_sat, m_sat)
        best = nt.norm(xx[cept,2]-xs.transpose()[cept])
        prev = best
        for t in range(13,17):
            dvv = np.zeros(dv.shape)
            d = dv[t_launch_indx]
            dvv[t_launch_indx] = dv[t_launch_indx]
            boost = 0.1
            print(t)
            for i in range(maxiter):
                xs, vs = ot.n_body_custom(mass, full_time, xx, vv, 1, dvv, False, x0_sat, v0_sat, m_sat)
                dist = nt.norm(xx[cept,2]-xs.transpose()[cept])#min(nt.norm(xx[indx,2] - xs.transpose(), ax = 1))
                if dist < best:
                    best = dist
                elif dist > prev:
                   boost = -0.5*boost
                print(dvv[t_launch_indx])
                d = d + boost
                dvv[t_launch_indx] = nt.rotate(d, -(2*np.pi)/360*0.25*t)
                print(dvv[t_launch_indx])
                print('Boost:',boost)
                print('dv:', dvv[t_launch_indx])
                print('TIME', full_time[t_launch_indx])
                prev = dist
                print('Closest approach:', dist)

        print('Best: ', best)
        return xs, dvv
        '''
    xss, dv = find_launch_sequence()
    a = nt.norm(xx[cept,2]-xss[cept - t_launch_indx])

    def inject():
        dv[cept] = -2.5*nt.unit_vector(xss.transpose()[cept])
        xs, vs = ot.n_body_custom(mass, full_time, xx, vv, 1, dv, False, x0_sat, v0_sat, m_sat)
        return xs
    xss = xss.transpose()
    for i in range(8):
        plt.plot(xx[:,i,0], xx[:, i, 1])
        plt.axis('equal')
    plt.plot(xss[0], xss[1], c = 'k')
    plt.scatter(xss[0, cept - t_launch_indx], xss[1, cept - t_launch_indx])
    plt.scatter(xx[cept, 2, 0], xx[cept, 2, 1])
    plt.show()

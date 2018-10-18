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

    steps = 20000
    t1 = 0.59
    tol = 0.0001
    time = np.linspace(0, t1, steps)
    import matplotlib.pyplot as plt

    # Generate super-orbits, hopefully only once.
    print('Simulating orbits.')
    mass = np.array([body.mass for body in sol.bodies])
    x0 = np.array([body.position for body in sol.bodies])
    v0 = np.array([body.velocity for body in sol.bodies])
    #xx, vv, nada, zipp = sol.find_orbits(time, 'sol', mass, x0, v0)
    '''
    np.save('saved/saved_params/xx_1mill.npy', xx)
    np.save('saved/saved_params/vv_1mill.npy', vv)
    '''
    print('Done')

    xx = np.load('saved/saved_params/xx_1mill.npy')
    vv = np.load('saved/saved_params/vv_1mill.npy')
    xx = xx.transpose()
    vv = vv.transpose()
    m_t = np.append(mass, m_sat)
    try:
        dw, tw, t_cept = ot.trajectory(m_t, xx, vv, steps, 1, -1, 2, 0, time, False, tol, 10000)
        print('Launch Parameters:', dw[0], tw[0], t_cept)
        #tw[0] = 0.42444244424442445; dw[0] = 2.2873197052846757
    except IndexError:
        print('No launch window found')

    # Simulating launch region with higher res.
    t_launch = tw[0] #0.4244386244386244/dt
    steps2 = 20000
    launc_duration = 0.001

    t1 = t_launch - t_launch/steps
    t2 = t_launch + 3*launc_duration
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
    dv = np.zeros(len(full_time))
    t_launch_indx = len(phase1) + 1
    cept = np.argwhere(abs(t_cept - full_time) < 0.001)[0,0]
    dv[t_launch_indx] = dw[0] # 2.25
    tol_dist = 0.0005

    print('Go for launch')

    def bazbaz(maxiter = 20):
        boost = 0.0001#dv[t_launch_indx]/100
        xs, vs = ot.n_body_custom(mass, full_time, xx, vv, 1, dv, False, x0_sat, v0_sat, m_sat)
        best = nt.norm(xx[cept,2]-xs.transpose()[cept])#min(nt.norm(xx[indx,2] - xs.transpose()[indx], ax = 1))
        prev = best
        for t in range(5):
            dvv = np.zeros(len(dv))
            print(t)
            for i in range(maxiter):
                dvv[t_launch_indx + t*500 + 500*15] = dv[t_launch_indx] + boost
                print('dv:', dv[t_launch_indx])
                print('TIME', full_time[t_launch_indx + t*2000])
                xs, vs = ot.n_body_custom(mass, full_time, xx, vv, 1, dvv, False, x0_sat, v0_sat, m_sat)
                dist = nt.norm(xx[cept,2]-xs.transpose()[cept])#min(nt.norm(xx[indx,2] - xs.transpose(), ax = 1))
                if dist < best:
                    best = dist
                elif dist < tol_dist:
                    break
                elif dist >= prev:
                   boost = -0.5*boost

                #dv[t_launch_indx] = dv[t_launch_indx] + boost
                prev = dist
                print('Closest approach:', dist)
        print('Best: ', best)
        return xs

    xss = bazbaz(5)
    a = nt.norm(xx[cept,2]-xss.transpose()[cept])#nt.norm(xx[indx,2] - xss.transpose()[indx], ax = 1)
    def inject():
        dv[cept] = 0.5
        xs, vs = ot.n_body_custom(mass, full_time, xx, vv, 1, dv, False, x0_sat, v0_sat, m_sat)
        return xs
    print(nt.norm(xss.transpose()[cept] - xss.transpose()[cept-1]))
    print(nt.norm(xx[cept,2]-xx[cept-1, 2]))
    xss = inject()

    xx = xx.transpose()
    for i in range(8):
        plt.plot(xx[0,i], xx[1,i])
        plt.axis('equal')
    plt.plot(xss[0], xss[1], c = 'k')
    plt.scatter(xss[0, cept], xss[1, cept])
    plt.scatter(xx[0,2, cept], xx[1,2, cept])
    plt.show()

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

    def find_orbits(self, time, ref_frame):
        steps = len(time)
        mass = np.array([body.mass for body in self.bodies])
        x0 = np.array([body.position for body in self.bodies])
        v0 = np.array([body.velocity for body in self.bodies])
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

    steps = 100000
    t1 = 0.58
    tol = 0.00001
    time = np.linspace(0, t1, steps)
    import matplotlib.pyplot as plt
    xx, vv, nada, zipp = sol.find_orbits(time, 'cm')

    mass = np.array([body.mass for body in sol.bodies])
    xx = xx.transpose()
    vv = vv.transpose()
    x0_sat = np.array([x0[0], y0[0]]) + vars.radius[0]*1000/vars.AU_tall
    v0_sat = np.array([vx0[0], vy0[0]])
    sat = Satellite(m_sat, x0_sat, v0_sat, 'MatSat')
    sol.addSatellite(sat)

    m_t = np.append(mass, m_sat)

    #dw, tw = ot.trajectory(m_t, xx, vv, steps, 1, -1, 2, 0, time, False, tol, 10000)
    #print(dw, tw)
    dv = np.zeros(len(time))
    dt = time[1]-time[0]
    t_launch = int(0.4244386244386244/dt)#tw[0]/dt)
    dv[t_launch] = 2.285416617717411#dw[0]

    xs, vs= ot.n_body_custom(mass, time, xx, vv, 1, dv, False, x0_sat, v0_sat, m_sat)

    while min(nt.norm(xx.transpose()[:,2] - xs)) > 0.001:
        dv[t_launch] = dv[t_launch] + 0.001
        xs, vs= ot.n_body_custom(mass, time, xx, vv, 1, dv, False, x0_sat, v0_sat, m_sat)
        print(min(nt.norm(xx.transpose()[:,2] - xs)))

    a = nt.norm(xx.transpose()[:,2] - xs)
    cept = np.unravel_index(np.argmin(a, axis=None), a.shape)

    xx = xx.transpose()

    for i in range(8):
        plt.plot(xx[0,i], xx[1,i])
        plt.axis('equal')
    plt.plot(xs[0], xs[1], '-.k')
    plt.scatter(xs[0, cept[0]], xs[1, cept[0]])
    plt.scatter(xx[0,2, cept[0]], xx[1,2, cept[0]])

    plt.axis('equal')
    plt.show()

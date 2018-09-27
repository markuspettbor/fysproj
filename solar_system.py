import numpy as np
import numtools as nt
import orbit_tools as ot

class SolarSystem:
    def __init__(self, num_planets, num_suns = 1):
        self.num_planets = num_planets
        self.planets = []
        self.suns = []
        self.bodies = []

    def addPlanet(self, planet):
        self.planets.append(planet)
        self.bodies.append(planet)

    def addSun(self, sun):
        self.suns.append(sun)
        self.bodies.append(sun)

    def find_orbits(self, t0, t1, steps, ref_frame):
        time = np.linspace(t0, t1, steps)
        mass = np.array([body.mass for body in self.bodies])
        x0 = np.array([body.position for body in self.bodies])
        v0 = np.array([body.velocity for body in self.bodies])
        names = [body.name for body in self.bodies]
        return ot.n_body_setup(mass, time, steps, x0, v0, ref_frame)

    def update_orbits(self, dt):
        pass
        bodies = self.suns + self.planets
        for i in range(self.num_planets):
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

if __name__ == '__main__':
    import variables
    m_star = variables.m_star
    m = variables.m
    x0 = variables.x0
    y0 = variables.y0
    vx0 = variables.vx0
    vy0 = variables.vy0
    radius = variables.radius
    names = ['Dum','og', 'Deilig', 'Juba', 'juba', 'Pizzatryne', 'Verdens ende']
    sol = SolarSystem(6, 1)
    sun = Sun(m_star, 0.001, np.array([0,0]), np.array([0,0]), 1000, 'Sol')
    sol.addSun(sun)


    for name, xx0, yy0, vvx0, vvy0, mass, r in zip(names, x0, y0, vx0, vy0, m, radius):
        x = np.array([xx0, yy0])
        v = np.array([vvx0, vvy0])
        planet = Planet(mass, r, x, v, name)
        sol.addPlanet(planet)

    xx, vv, nada, zipp = sol.find_orbits(0, 10, 1000, 'not_cm')

    import matplotlib.pyplot as plt
    for i in range(6):
        plt.plot(xx[0,i], xx[1,i])
    plt.axis('equal')
    plt.show()

    #def potential_energy(system_masses, system_x, target_x, steps, body_index = 0):
    mass = np.array([body.mass for body in sol.bodies])
    xx = xx.transpose()
    ot.trajectory(mass, xx, 1000, 2, 1, 3, 0 )

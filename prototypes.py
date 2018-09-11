import numpy as np
from solar_system import SolarSystem, Planet, Sun
import pygame as pg
import variables as vars
import numtools as nt

class Screen(object):
    def __init__(self, w, h):
        self.w = w
        self.h = h
        pg.init()
        pg.display.set_mode((self.w, self.h))

    def update(self):
        pass

    def close(self):
        pg.display.quit()


def gravity(m1, m2, x):
    return -G*m1*m2/nt.norm(x)**3*x

if __name__ == '__main__':
    width = 1920 # Can't remember resolution things
    height = 1280
    test = Screen(width, height)
    testrun = True
    surf = pg.display.get_surface()
    user = vars.seed
    solar_system = vars.solar_system
    n = vars.n
    x0 = vars.x0
    y0 = vars.y0
    vx0 = vars.vx0
    vy0 = vars.vy0
    a = vars.a
    e = vars.e
    theta0 = vars.theta0
    psi0 = vars.psi0
    radius = vars.radius
    m_star = vars.m_star
    m = vars.m
    G = vars.G

    planet_names = ['Dum','og', 'Deilig', 'Juba', 'Juba', 'Kampala', 'Pizzatryne', 'Verdens ende']
    sun_name = 'pleb'
    fullmugg = SolarSystem(n, 1)
    sunsun = Sun(m_star, 0.001, np.array([0,0]), 100000, 'sunsun')
    fullmugg.addSun(sunsun)

    # Init solar system
    for name, x0, y0, vx0, vy0, m, r in zip(planet_names, x0, y0, vx0, vy0, m, radius):
        x = np.array([x0, y0])
        v = np.array([vx0, vy0])/2
        planet = Planet(m, r, x, v, name)
        fullmugg.addPlanet(planet)

    dt = 0.0001
    x = np.array([planet.position for planet in fullmugg.planets])
    v = np.array([planet.velocity for planet in fullmugg.planets])
    count = 0
    testrun = True
    scale = 100
    center = 700
    surf.fill((255,255,255))
    while testrun:
        pg.draw.circle(surf, (0,0,0), (center, center), 10)
        for planet, index in zip(fullmugg.planets, range(n)):
            m = planet.mass
            acceleration = lambda r: gravity(m_star, m, r)/m
            x[index], v[index] = nt.leapfrog_simple(x[index], v[index], dt, acceleration)
        if count % 10 == 0:
            for xx, yy in zip(x[:, 0], x[:,1]):
                pg.draw.circle(surf, (0,0,0), (int(xx*scale) + center, int(yy*scale)+center), 0)
            pg.display.flip()
            #surf.fill((255,255,255))

        for event in pg.event.get():
            if event.type == pg.KEYDOWN and event.key == pg.K_q:
                test.close()
                testrun = False
        count +=1

'''
            def acc(r):
                mask = np.arange(n) != index
                return np.sum(gravity(m_star, m[mask], r[index] - r[mask])/m[index])

'''

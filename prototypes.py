import numpy as np
import matplotlib.pyplot as plt
from solar_system import SolarSystem, Planet, Sun
import pygame as pg
import variables as vars
import numtools as nt
from numba import jit

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

    names = ['Matsat','sunsun', 'Dum','og', 'Deilig', 'Juba', 'juba', 'Pizzatryne', 'Verdens ende']
    sun_name = 'pleb'

    fullmugg = SolarSystem(n+1, 1)
    sunsun = Sun(m_star, 0.001, np.array([0,0]), 100000, 'sunsun')
    fullmugg.addSun(sunsun)

    masses = np.concatenate((np.array([m_star]), m))

    radius = np.concatenate((np.array([vars.radius_star]), radius))
    x0 = np.concatenate((np.array([0]), x0))
    y0 = np.concatenate((np.array([0]), y0))
    vx0 = np.concatenate((np.array([0]), vx0))
    vy0 = np.concatenate((np.array([0]), vy0))

    initial_rocket_mass = 4.097e+10/1.989e30
    m_sat = 1100/1.989e30

    m_rock = np.array([(initial_rocket_mass + m_sat)])
    x0_sat = np.array([x0[1]])
    y0_sat = np.array([y0[1]+ radius[1]*1000/vars.AU_tall])
    vx0_sat = np.array([vx0[1]])
    vy0_sat = np.array([vy0[1]]) + 11000*60*60*24*365/vars.AU_tall
    radius_sat = np.array([2/vars.AU_tall])

    masses = np.concatenate((m_rock, masses))
    x0 = np.concatenate((x0_sat, x0))
    y0 = np.concatenate((y0_sat, y0))
    vx0 = np.concatenate((vx0_sat, vx0))
    vy0 = np.concatenate((vy0_sat, vy0))
    radius = np.concatenate((radius_sat, radius))

    # Init solar system
    for name, xx0, yy0, vvx0, vvy0, mass, r in zip(names, x0, y0, vx0, vy0, masses, radius):
        x = np.array([xx0, yy0])
        v = np.array([vvx0, vvy0])
        planet = Planet(mass, r, x, v, name)
        fullmugg.addPlanet(planet)

    dt = 0.0001
    count = 0
    testrun = True
    scale = 300
    center = 700

    width = 1920 # Can't remember resolution things
    height = 1280
    test = Screen(width, height)
    testrun = True
    surf = pg.display.get_surface()

    surf.fill((255,255,255))
    font = pg.font.SysFont("comicsansms", 15)

    x00 = np.array([x0, y0])
    v00 = np.array([vx0, vy0])
    x0 = x00.transpose()
    v0 = v00.transpose()
    n = len(masses)
    steps = 1
    cm  = np.zeros((steps, 2))
    vcm = np.zeros((steps, 2))
    xx = np.zeros((1, n, 2))
    vv = np.zeros((1, n, 2))
    xx[0] = x0
    vv[0] = v0
    x = x0
    v = v0

    import orbit_tools as ot

    # All params collected from test_tools and launch.py
    indx = 0

    engine_boxes = 3.909e+15
    fuel_consumed = 5.83160e-14/(1.989e30)*engine_boxes*365*24*60*60  # solmasses/yr
    dv = 11000/vars.AU_tall/(5.630*60)*(60*60*24*365) # Approximate delta-v AU/yr**2
    print(fuel_consumed*1.989e30/(365*60*60*24), dv*vars.AU_tall/(60*60*24*365))
    while testrun:
        #pg.draw.circle(surf, (0,0,0), (center, center), 2)
        for k in range(1):
            x = np.copy(xx[k])
            v = np.copy(vv[k])
            for i in range(n):
                acc = lambda r: ot.system_acceleration(masses, r, i, n)
                x[i] = xx[k,i] + vv[k,i]*dt + 0.5*acc(xx[k])*dt**2
            for i in range(n):
                acc = lambda r: ot.system_acceleration(masses, r, i, n)
                v[i] = vv[k,i] + 0.5*(acc(xx[k])+ acc(x))*dt

            cm = ot.center_of_mass(masses, x)
            xx[k] = x
            vv[k] = v
        if count % 1 == 0:
            for xd, yd , planet in zip(x[:,0], x[:,1], fullmugg.planets):
                #pg.draw.circle(surf, (0,0,0), (int(0 + center) , int(0 + center)), 3)
                xd = xd - cm[0] - xx[0, indx, 0]
                yd = yd - cm[1] - xx[0, indx, 1]

                if planet.name == 'Matsat':
                    col = (100, 190, 100)
                else:
                    col = (0, 0, 0)
                pg.draw.circle(surf, col, (int(xd*scale) + center, int(yd*scale)+center), 3)
                strstr = planet.name + ' ' + '%.2e' %(planet.mass)

                text = font.render(strstr, True, col)
                surf.blit(text, (int(xd*scale) + center, int(yd*scale)+center))
            pg.display.flip()
        surf.fill((255,255,255))

        for event in pg.event.get():
            if event.type == pg.KEYDOWN:
                if event.key == pg.K_q:
                    test.close()
                    testrun = False
                if event.key == pg.K_b:
                    dt = dt/2
                if event.key == pg.K_n:
                    dt = dt*2
                if event.key == pg.K_w:
                    scale = scale*2
                    surf.fill((255,255,255))
                if  event.key == pg.K_s:
                    scale = scale/2
                    surf.fill((255,255,255))
                if event.key == pg.K_0:
                    indx = 0
                if event.key == pg.K_1:
                    indx = 1
                if event.key == pg.K_2:
                    indx = 2
                if event.key == pg.K_3:
                    indx = 3
                if event.key == pg.K_4:
                    indx = 4
                if event.key == pg.K_5:
                    indx = 5
                if event.key == pg.K_6:
                    indx = 6
                if event.key == pg.K_7:
                    indx = 7

        keys = pg.key.get_pressed()
        if fullmugg.planets[0].mass <= m_sat:
            #dv = 0
            #fuel_consumed = 0
            #fullmugg.planets[0].mass = m_sat
            masses[0] = m_sat

        if keys[pg.K_i]:
            vv[0,0] = vv[0,0]+ np.array([0, -dv])
            print(fuel_consumed*dt, dt, fullmugg.planets[0].mass)
            fullmugg.planets[0].mass -= fuel_consumed*dt
            masses[0] -= fuel_consumed*dt

        if keys[pg.K_k]:
            vv[0,0] = vv[0,0]+ np.array([0, dv])
            fullmugg.planets[0].mass -= fuel_consumed*dt
            masses[0] -= fuel_consumed*dt

        if keys[pg.K_l]:
            vv[0,0] = vv[0,0]+ np.array([dv, 0])
            fullmugg.planets[0].mass -= fuel_consumed*dt
            masses[0] -= fuel_consumed*dt

        if keys[pg.K_j]:
            vv[0,0] = vv[0,0] + np.array([-dv, 0])
            fullmugg.planets[0].mass -= fuel_consumed*dt
            masses[0] -= fuel_consumed*dt

        if keys[pg.K_m]:
            fullmugg.planets[0].mass +=  0.0001
            masses[0] += 0.0001

        count +=1

    def acc(r):
        mask = np.arange(n) != index
        return np.sum(gravity(m_star, m[mask], r[index] - r[mask])/m[index])

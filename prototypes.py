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
    sunsun = Sun(m_star, 0.001, np.array([0,0]), np.array([0,0]), 100000, 'sunsun')
    fullmugg.addSun(sunsun)

    masses = np.concatenate((np.array([m_star]), m))

    radius = np.concatenate((np.array([vars.radius_star]), radius))
    x0 = np.concatenate((np.array([0]), x0))
    y0 = np.concatenate((np.array([0]), y0))
    vx0 = np.concatenate((np.array([0]), vx0))
    vy0 = np.concatenate((np.array([0]), vy0))

    initial_rocket_mass = 1.456e+05/1.989e30
    m_sat = 1100/1.989e30

    m_rock = np.array([m_sat])
    x0_sat = np.array([x0[1]])
    y0_sat = np.array([y0[1]+ radius[1]*1000/vars.AU_tall])
    vx0_sat = np.array([vx0[1]])
    vy0_sat = np.array([vy0[1]]) #+ 10500*60*60*24*365/vars.AU_tall
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

    dt = 0.00001
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
    xx = x0
    vv = v0
    x = x0
    v = v0

    import orbit_tools as ot

    # All params collected from test_tools and launch.py
    indx = 0
    frameno = 1

    engine_boxes = 1.368e+17
    fuel_consumed = 5.83160e-14/(1.989e30)*engine_boxes*365*24*60*60  # solmasses/yr
    force = 2.55846e-10*engine_boxes/1.989e30/vars.AU_tall*(365*24*60*60)**2
    dv = force/masses[0]*dt

    t_boost = 0.42193021930219304
    tcount = 0
    dv_boost =  2.28271470381530
    not_launched = True
    closest = 10
    while testrun:
        #pg.draw.circle(surf, (0,0,0), (center, center), 2)

        x = np.copy(xx)
        v = np.copy(vv)
        for i in range(n):
            acc = lambda r: ot.system_acceleration(masses, r, i, n)
            x[i] = xx[i] + vv[i]*dt + 0.5*acc(xx)*dt**2
        for i in range(n):
            acc = lambda r: ot.system_acceleration(masses, r, i, n)
            v[i] = vv[i] + 0.5*(acc(xx)+ acc(x))*dt

        cm = ot.center_of_mass(masses, x)
        xx = x
        vv = v

        tcount += dt
        if tcount > t_boost and not_launched:
            dt = dt/10000000
            xx[0] = np.copy(xx[2]) + vars.radius[0]*1000/vars.AU_tall*nt.unit_vector(vv[2])
            vv[0] = np.copy(vv[2]) + nt.unit_vector(vv[2])*dv_boost
            not_launched = False
            print('LAUNCHED!')

        elif not_launched:
            xx[0] = np.copy(xx[2]) + vars.radius[0]*1000/vars.AU_tall*nt.unit_vector([vv[2]])
            vv[0] = np.copy(vv[2])


        small = nt.norm(xx[0]- xx[3])
        if small < closest:
            closest = small
            print('Closest Approach: ', closest)


        if count % 1 == 0:
            for xd, yd , planet in zip(x[:,0], x[:,1], fullmugg.planets):
                #pg.draw.circle(surf, (0,0,0), (int(0 + center) , int(0 + center)), 3)
                xd = xd - cm[0] - xx[indx, 0]*frameno
                yd = yd - cm[1] - xx[indx, 1]*frameno

                if planet.name == 'Matsat':
                    col = (100, 190, 100)
                else:
                    col = (0, 0, 0)
                pg.draw.circle(surf, col, (int(xd*scale) + center, int(yd*scale)+center), 0)
                #strstr = planet.name + ' ' + '%.2e' %(planet.mass)

                #text = font.render(strstr, True, col)
                #surf.blit(text, (int(xd*scale) + center, int(yd*scale)+center))
            pg.display.flip()
        #surf.fill((255,255,255))

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
                    frameno = 1
                    surf.fill((255,255,255))
                if event.key == pg.K_1:
                    surf.fill((255,255,255))
                    indx = 1
                    frameno = 1
                if event.key == pg.K_2:
                    surf.fill((255,255,255))
                    indx = 2
                    frameno = 1
                if event.key == pg.K_3:
                    surf.fill((255,255,255))
                    indx = 3
                    frameno = 1
                if event.key == pg.K_4:
                    surf.fill((255,255,255))
                    indx = 4
                    frameno = 1
                if event.key == pg.K_5:
                    surf.fill((255,255,255))
                    indx = 5
                    frameno = 1
                if event.key == pg.K_6:
                    surf.fill((255,255,255))
                    indx = 6
                    frameno = 1
                if event.key == pg.K_7:
                    surf.fill((255,255,255))
                    indx = 7
                    frameno = 1
                if event.key == pg.K_9:
                    surf.fill((255,255,255))
                    frameno = 0

        keys = pg.key.get_pressed()
        dv = force/masses[0]*dt
        dv = 0.01
        '''
        if masses[0] - m_sat < fuel_consumed*dt and masses[0] - m_sat > 0:
            boost_time = (masses[0]-m_sat)/fuel_consumed
            dv = force/masses[0]*boost_time

        if fullmugg.planets[0].mass <= m_sat:
            dv = 0
            fuel_consumed = 0
            fullmugg.planets[0].mass = m_sat
            masses[0] = m_sat
        '''
        if keys[pg.K_i]:
            vv[0,0] = vv[0,0]+ np.array([0, -dv])
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

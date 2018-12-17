import numpy as np
import prob_dist as pd
import matplotlib.pyplot as plt
import numtools as nt
import matplotlib
import physicals as ph
import pygame as pg
import random
import rocket as rock
from launch import launch
import variables as vars
# This is our code

class Screen(object):
    # Pygame needed for visualization.
    def __init__(self, w, h):
        self.w = w
        self.h = h
        pg.init()
        pg.display.set_mode((self.w, self.h))

    def update(self):
        pass

    def close(self):
        pg.display.quit()

def tasks():
    # Task 1
    np.random.seed(69)
    n = 10**5
    x0 = np.linspace(-2.5, 2.5, 10000)*10**4
    boltzmax = pd.BoltzmannMaxwell(T = 1e4, N = n) # Maxwell Boltzmann distribution
    plt.plot(x0, boltzmax(x0), '-k', linewidth = 0.8)
    plt.xlabel('$v_x$ [m/s]', size = 14)
    plt.ylabel('$P(v_x)$', size  = 14)
    plt.title('Velocity Distribution, $P(v_x)$')
    plt.show()

    # Task 2
    x1 = np.linspace(5, 30, 10000)*10**3
    area = nt.integrate(lambda x: boltzmax(x1), x1, teacher_is_strict = False)
    print('Area:', area, n*area)

    # Task 3
    x2 = np.linspace(0, 3, 10000)*10**4
    absolute_boltzmax = pd.AbsoluteBoltzmannMaxwell(T = 1e4, N = n)
    plt.plot(x2, absolute_boltzmax'''
k =  1.380 648 52 x 10-23 J K-1
https://physics.nist.gov/cgi-bin/cuu/Value?k|search_for=boltzmann

Na = 6.022 140 857 x 1023 mol-1
https://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=avogadro

mass H_2 = 2.016 g/mol
https://pubchem.ncbi.nlm.nih.gov/compound/Hydrogen
'''.absolute_density(x2), '-k', linewidth = 0.8)
    plt.xlabel('$v$ [m/s]', size = 14)
    plt.ylabel('$P(v)$', size = 14)
    plt.title('Absolute Velocity Distribution, $P(v)$')
    plt.show()

def chamber():
    width = 800 # Resolution
    height = 800
    test = Screen(width, height)
    testrun = True
    surf = pg.display.get_surface()
    box_width = 1e-6
    x0 = 0
    y0 = box_width
    z0 = box_width # Box dimensions
    n_part = 5000
    temp = 3000
    h2 = ph.Gas(n_part, temp, -box_width/2, box_width/2) # H2 gas
    p = h2.position.transpose() + np.array([0, y0, z0])
    p = p.transpose()
    v = h2.velocity
    bigbox = ph.Box(x0, y0, z0 , box_width, hole_index = 4, hole_width = box_width/2) # creates box
    walls = bigbox.walls # Create walls objects, that can be collided with
    dt = 1e-12 #må stemme overens med dimensjoner på box og hastigheter
    white = (255, 255, 255, 255)
    count = 0
    force = 0
    fuel_used = 0
    tot_force = 0
    engine = rock.Engine(1000,100)
    intervals = 1e3
    scale = 400

    while testrun:
        for wall in walls:
            # Determine collisions
            v = wall.check_collision(p, v)
            force_temp, fuel_used_temp = engine.particles_escaped(wall.escaped_particles, \
                                             wall.escaped_velocity, dt)
            force += force_temp
            fuel_used += fuel_used_temp
            wall.reset_escaped_particles()

        v, p = nt.euler_cromer_simple(p, v, dt)
        if count % 20 == 0:
            # Visualize chamber
            for x, y in zip(p[1], p[2]):
                pg.draw.circle(surf, white, (int(x*scale/box_width), int(y*scale/box_width)), 0)
            pg.display.flip()
            surf.fill((0,0,0))
        count += 1
        tot_force += force
        force = 0
        if count % int(intervals) == 0:
            # Determine engine/ launch params
            print('Force %19.5e [N]' % (tot_force/intervals))
            print('Fuel consumed %11.5e [kg/s]' % (fuel_used/intervals/dt))
            force_box = tot_force/intervals
            fuel_consumed_box_per_sec = fuel_used/intervals/dt
            tot_force = 0
            fuel_used = 0
            testrun = False
        for event in pg.event.get():
            if event.type == pg.KEYDOWN and event.key == pg.K_q:
                test.close()
                testrun = False
    area = (box_width/2)**2
    simulated_pressure = force_box*2/area #force * 2 is force delivered at collision...
    ideal_pressure = n_part * vars.k * temp / box_width**3
    compare = simulated_pressure/ideal_pressure
    print('sim/comp =', compare)
    def save_values():
        np.save('saved/engine/force_box.npy', force_box)
        np.save('saved/engine/fuel_consumed_box_per_sec.npy', fuel_consumed_box_per_sec)
    #save_values()
    #launch(force_box, fuel_consumed_box_per_sec, testing = True)
if __name__ == '__main__':
    tasks()
    chamber()

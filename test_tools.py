import prob_dist as pd
import physicals as ph
import numtools as nt
import pygame as pg
import random
import numpy as np
import rocket_parts as rp
import prototypes
from launch import launch

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

if __name__ == '__main__':
    width = 800 # Can't remember resolution things
    height = 800
    test = Screen(width, height)
    testrun = True
    surf = pg.display.get_surface()

    box_width = 1e-6
    x0 = 0
    y0 = box_width#int(width/2)
    z0 = box_width#int(height/2)

    h2 = ph.Gas(50000, 3000, -box_width/2, box_width/2)
    p = h2.position.transpose() + np.array([0, y0, z0])
    p = p.transpose()
    v = h2.velocity

    # Alternative box:


    bigbox = prototypes.Box(x0, y0, z0 , box_width, hole_index = 4, hole_width = box_width/2)
    walls = bigbox.walls

    dt = 1e-12 #må stemme overens med dimensjoner på box og hastigheter
    white = (255, 255, 255, 255)
    count = 0
    force = 0
    fuel_used = 0
    tot_force = 0
    engine = rp.Engine(1000,100)
    intervals = 1e3
    scale = 400


    while testrun:
        for wall in walls:
            #pg.draw.circle(surf, white, (int(wall.center[1]), int(wall.center[2])), 2)
            v = wall.check_collision(p, v)
            force_temp, fuel_used_temp = engine.particles_escaped(wall.escaped_particles, \
                                             wall.escaped_velocity, dt)
            force += force_temp; fuel_used += fuel_used_temp
            wall.reset_escaped_particles()

        v, p = nt.euler_cromer_simple(p, v, dt)
        if count % 20 == 0:
            for x, y in zip(p[1], p[2]):
                pg.draw.circle(surf, white, (int(x*scale/box_width), int(y*scale/box_width)), 0)
            pg.display.flip()
            surf.fill((0,0,0))
        count += 1
        #print(force)
        tot_force += force
        force = 0
        if count % int(intervals) == 0:
            print('Force %19.5e [N]' % (tot_force/intervals))
            print('Fuel consumed %11.5e [1/s]' % (fuel_used/intervals/dt))
            force_box = tot_force/intervals
            fuel_consumed_box_per_sec = fuel_used/intervals/dt
            tot_force = 0
            fuel_used = 0
            testrun = False
        for event in pg.event.get():
            if event.type == pg.KEYDOWN and event.key == pg.K_q:
                test.close()
                testrun = False
    launch(force_box, fuel_consumed_box_per_sec, testing = True)

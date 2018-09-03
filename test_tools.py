import prob_dist as pd
import physicals as ph
import numtools as nt
import pygame as pg
import random
import numpy as np

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
    height = 600
    test = Screen(width, height)
    testrun = True
    surf = pg.display.get_surface()

    x0 = 0
    y0 = int(width/2)
    z0 = int(height/2)
    box_width = 400

    h2 = ph.Gas(10, 3000, -box_width/2, box_width/2)
    p = h2.position.transpose() + np.array([0, y0, z0])
    p = p.transpose()
    v = h2.velocity

    # Alternative box:
    import prototypes

    bigbox = prototypes.Box(x0, y0, z0 , box_width, hole_index = 4)
    walls = bigbox.walls

    dt = 0.00005
    white = (255, 255, 255, 255)
    count = 0

    while testrun:

        for wall in walls:
            pg.draw.circle(surf, white, (int(wall.center[1]), int(wall.center[2])), 2)
            v = wall.check_collision(p, v)
        v, p = nt.euler_cromer_simple(p, v, dt)
        if count % 1 == 0:
            for x, y in zip(p[1], p[2]):
                pg.draw.circle(surf, white, (int(x), int(y)), 0)
            pg.display.flip()
            surf.fill((0,0,0))
        count += 1
        for event in pg.event.get():
            if event.type == pg.KEYDOWN and event.key == pg.K_q:
                test.close()
                testrun = False

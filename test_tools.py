import prob_dist as pd
import physicals
import numtools as nt
import pygame as pg
import random

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
    width = 1920 # Can't remember resolution things
    height = 1080
    test = Screen(width, height)
    testrun = True
    surf = pg.display.get_surface()
    h2 = physicals.Gas(100000, 1000)
    px = h2.position[:, 0] + 500
    py = h2.position[:, 1] + 500
    vx = h2.velocity[:, 0]
    vy = h2.velocity[:, 1]
    dt = 0.001
    white = (255, 255, 255, 255)
    count = 0
    while testrun:
        pg.draw.circle(surf, (255, 255, 255, 255), (1000, 1000), 10)
        vx, px = nt.euler_cromer_simple(px, vx, dt)
        vy, py = nt.euler_cromer_simple(py, vy, dt)
        for x, y in zip(px, py):
            pg.draw.circle(surf, (255, 255, 255, 255), (int(x), int(y)), 1)
        if count % 1 == 0:
            pg.display.flip()
            surf.fill((0,0,0))

        count += 1

        for event in pg.event.get():
            if event.type == pg.KEYDOWN and event.key == pg.K_q:
                test.close()
                testrun = False

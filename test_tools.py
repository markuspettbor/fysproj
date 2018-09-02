import prob_dist as pd
import physicals
import numtools as nt
import pygame as pg
import random
import numpy as np
from prototypes import Wall

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
    h2 = physicals.Gas(1000, 1000)
    px = h2.position[:, 0] + 500
    py = h2.position[:, 1] + 500
    vx = h2.velocity[:, 0]
    vy = h2.velocity[:, 1]

    p = h2.position.transpose() + 500
    v = h2.velocity.transpose()

    n = np.array([0, 1, 0])
    axis = 0
    sign = 1
    center = np.array([700, 500 , 0])
    w1 = Wall(n, axis, sign, center)

    n = np.array([0, -1, 0])
    axis = 0
    sign = -1
    center = np.array([300, 500 , 0])
    w2 = Wall(n, axis, sign, center)

    n = np.array([1, 0, 0])
    axis = 1
    sign = 1
    center = np.array([500, 700 , 0])
    w3 = Wall(n, axis, sign, center)

    n = np.array([-1, 0, 0])
    axis = 1
    sign = -1
    center = np.array([500, 300 , 0])
    w4 = Wall(n, axis, sign, center)

    walls = [w1, w2, w3, w4]

    dt = 0.0001
    white = (255, 255, 255, 255)
    count = 0

    while testrun:
        pg.draw.circle(surf, white, (700, 500), 10)
        pg.draw.circle(surf, white, (300, 500), 10)
        pg.draw.circle(surf, white, (500, 700), 10)
        pg.draw.circle(surf, white, (500, 300), 10)
        for wall in walls:
            v = wall.check_collision(p, v)
        v, p = nt.euler_cromer_simple(p, v, dt)
        for x, y in zip(p[0], p[1]):
            pg.draw.circle(surf, (255, 255, 255, 255), (int(x), int(y)), 1)
        if count % 1 == 0:
            pg.display.flip()
            surf.fill((0,0,0))

        count += 1

        for event in pg.event.get():
            if event.type == pg.KEYDOWN and event.key == pg.K_q:
                test.close()
                testrun = False

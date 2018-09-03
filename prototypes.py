import numpy as np


def test():
    n = np.array([0, 0, 1])
    axis = 2
    sign = 1
    center = np.array([0, 0 , 5])
    num = 100
    w1 = Wall(n, axis, sign, center, hole_width = 100)
    pos = np.random.randint(0, 10, size = (3, num))
    vel = np.ones((3, num))
    w1.boundary(pos)
test()

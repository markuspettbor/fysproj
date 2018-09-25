import numpy as np

class Engine:
    def __init__(self, length, width, height):
        pass

    def make_wall(self):
        pass

    def detect_collision(self):
        pass

#[x, y, vx, vy, mass, radius, force]
import matplotlib.pyplot as plot
import variables as vars
t = np.load('saved_orbits/launch_resolution/time.npy')
print(t[0:10]*365*24*60)
print(max(t)*365*24*60)
print(t.shape)

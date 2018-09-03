import random
import physicals as ph

class Box:

    def __init__(self, center, dim, walls):
        '''
        Assumes center is tuple of coordinates (x,y,x) for box center.
        Assumes size is a tuple containing box size (width, height, depth)
        '''
        self.center = center
        self.dim = dim
        self.walls = walls
        self.exited_particles = 0

    def make_box(self):
        pass

import random
import physicals as ph

class Box:

    def __init__(self, center, dim):
        '''
        Assumes center is tuple of coordinates (x,y,x) for box center.
        Assumes size is a tuple containing box size (width, height, depth)
        '''
        self.center = center
        self.dim = dim
        self.w = self.dim[0]
        self.h = self.dim[1]
        self.d = self.dim[2]
        sel.walls = {}

    def make_box(self):
        sides = ('top', 'bottom', 'left', 'right', 'front', 'back')
        index = ((self.w/2, 1), (-self.w/2, 1), (self.h/2, 2), (-self.h/2, 2), (self.d/2, 0), (-self.d/2, 0))
        for wall, location in zip(sides, index):
            self.walls[wall] = location


def test_box():
    test = Box((0,0,0), (1, 2, 3))
    test.make_box()
    print(test.walls)

test_box()

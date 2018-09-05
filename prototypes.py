import numpy as np
import physicals

class Box:
    def __init__(self, dx, dy, dz, width, hole_index = 5, hole_width = 200):
        '''
        Creates the simplest of all boxes.
        '''
        normals = np.array([[1,0,0], [-1, 0, 0], [0, 1, 0],\
                            [0, -1, 0], [0, 0, 1], [0, 0, -1]])
        axes = np.array([0, 0, 1, 1, 2, 2])
        signs = np.array([1, -1, 1, -1, 1, -1])
        x1 = np.array([dx + width/2, dy, dz])
        x2 = np.array([dx - width/2, dy, dz])
        x3 = np.array([dx, dy + width/2, dz])
        x4 = np.array([dx, dy - width/2, dz])
        x5 = np.array([dx, dy, dz + width/2])
        x6 = np.array([dx, dy, dz - width/2])
        centers = [x1, x2, x3, x4, x5, x6]
        self.walls = []
        holes = np.zeros(6)
        holes[hole_index] = hole_width
        self.make_box(normals, axes, signs, centers, holes)

    def add_wall(self, normal, axis, sign, center, hole):
        return physicals.Wall(normal, axis, sign, center, hole_width = hole)

    def make_box(self, norms, axes, signs, centers, holes):
        for n, ax, sign, cntr, hole in zip(norms, axes, signs, centers, holes):
            self.walls.append(self.add_wall(n, ax, sign, cntr, hole))


def test():
    bigbox = Box(500, 500, 500, 200)
    for wall in bigbox.walls:
        print(wall.hole_width)

test()

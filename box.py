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

    def check_collision(self, position, velocity):
        for wall in self.walls:
            normal = wall.direction()
            velocity = np.where(wall.boundary(position),
                        velocity*(normal), velocity)
        return velocity

    def move_particles(self):
        pass

    def detect_exit(self, position):
        for wall in self.walls:
            try:
                hole = wall.hole(position)
                temp = np.where(position > hole)
                escaped = np.count_nonzero(temp)
                self.exited_particles += escaped
                position = np.where(position > hole, position, self.move_particles)
            except:
                print('Wall has no hole')
        return position


def test_box():
    pass

test_box():

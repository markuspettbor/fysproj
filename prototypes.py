import numpy as np

class Wall():
    def __init__(self, normal_vector, axis, sign, center, corners = None, hole_width = 0, molecule_moment = 0):
        self.normal_vector = normal_vector
        self.axis = axis
        self.sign = sign
        self.center = center
        self.corners = corners
        self.hole_width = hole_width
        self.unit_normal = normal_vector/np.linalg.norm(normal_vector)
        self.escaped = 0
        self.dp = 0
        self.p = molecule_moment

    def check_collision(self, position, velocity):
        velocity[self.axis] = np.where(self.boundary(position),\
                            -velocity[self.axis], velocity[self.axis])
        return velocity

    def boundary(self, position):
        '''
        Assumes position is a numpy array containing positions (x,y,z),
        with dimension (3, N).
        First checks if any particle is outside the wall along a given axis.
        If there is a hole, and particles are outside, it checks to see if those
        particles are within the bounds of the hole.

        BUGTOFIX: If the position is zero, the check can trigger in a rubbish way.
        '''
        outside = self.sign*position[self.axis] >= self.sign*self.center[self.axis]
        # Find out which particles are outside. Then, check among those that are outside, if they are in hole!
        match = np.any(outside)
        if self.hole_width is not 0 and match:
            mask_index = np.arange(3) != self.axis # retrieve relevant axes
            grid = np.where(outside, position[mask_index], False)
            in_hole = np.abs(self.center[mask_index] - grid.transpose()) <= self.hole_width/2
            in_hole = in_hole.transpose()
            esc = in_hole[0]*in_hole[1]*outside
            self.escaped += np.count_nonzero(esc)
        return outside

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
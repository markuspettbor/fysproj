import numpy as np
#import numba

def integrate(func, t, teacher_is_strict = True):
    '''
    Function for integration.
    Assumes func is a function that evaluates a numpy array of x-values.
    If teacher_is_strict is True, the midpoint method is used.
    If not, the built-in numpy trapz-function is used.
    Returns sum, which is integrated evaluated over the range of x.
    '''
    try:
        assert teacher_is_strict
        dt = np.diff(t)
        t_mid = t[:-1] + dt/2
        sum = np.sum(func(t_mid)*dt)
    except AssertionError:
        sum = np.trapz(func(t), t)
    return sum

def norm(vector, ax = 0):
    return np.linalg.norm(vector, axis = ax)

def unit_vector(vector, ax = 0):
    return norm(vector, ax)

def euler_cromer_simple(x, v, dt, acc = 0):
    v = v + dt*acc
    x = x + dt*v
    return v, x


def leapfrog(x0, v0, t, acc):
    '''
    Calculates integral using a Leapfrog method.
    x0, y0 are initial values, t is a vector of intervals. Note that it
    is assumed that dt = t(1) - t(0) is constant. It has to be, in order for
    energy to be conserved, apparently.
    acc is a method or function that evaluates the acceleration at a given time
    t.
    Assumes x0, v0 are of same length
    '''
    x = np.zeros((len(t), len(x0)))
    v = np.zeros((len(t), len(v0)))
    x[0] = x0
    v[0] = v0
    dt = t[1] - t[0]
    for i in range(len(t)-1):
        x[i+1] = x[i] + v[i]*dt + 0.5*acc(x[i], t[i])*dt**2
        v[i+1] = v[i] + 0.5*(acc(x[i], t[i]) + acc(x[i+1], t[i+1]))*dt
    return np.transpose(x), np.transpose(v)

def leapfrog_simple(x0, v0, dt, acc):
    x = x0 + v0*dt + 0.5*acc(x0)*dt**2
    v = v0 + 0.5*(acc(x0)+ acc(x))*dt
    return x, v

def euler_fuel_consumption(speed, mass, force, consumption, dt = 0.001):
    a = force/mass
    speed = speed + a*dt
    mass = mass - consumption*dt
    return speed, mass

def rotate(vector, angle):
    # Rotation matrix
    x1 = np.cos(angle)
    x2 = np.sin(angle)
    rotmatrix = np.array([[x1, -x2], [x2, x1]])
    return np.dot(vector, rotmatrix)


def angle_between(v1, v2):
    """ https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    Returns the angle in radians between vectors 'v1' and 'v2':: """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def create_radial_velocity(v, v_pec_r, i):
    vr = v[1]*np.sin(i) #ROW OR COL???
    v_real = vr + v_pec_r
    #noice = np.random.normal(0, max(vr)/5, len(vr)) #mu sigma length
    return noiceify(v_real, 0, max(vr)/5)

def least_squares(function, vars):
    pass
    return noiceify(v_real, 0, max(vr)/5)
    #noice = np.random.normal(0, max(vr)/5, len(vr)) #mu sigma length
    #return (v_real + noice)

def noiceify(x, mu, sig):
    return (x + np.random.normal(mu, sig, len(x))) #mu sigma length

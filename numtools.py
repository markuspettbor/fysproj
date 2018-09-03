import numpy as np

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

def angle_between(v1, v2):
    """ https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    Returns the angle in radians between vectors 'v1' and 'v2':: """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

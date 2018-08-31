import numpy as np
import itertools
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

def combos(vectors)

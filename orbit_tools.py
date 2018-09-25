import numpy as np
import numtools as nt
import variables as vars

def gravity(m1, m2, x):
    return -vars.G*m1*m2/nt.norm(x)**3*x

def kepler3(m1, m2, a):
    return np.sqrt(4*np.pi**2/(vars.G*(m1+m2))*a**3)

def orbit(x0, v0, acc, t):
    '''
    Assumes x0, v0, are arrays containing [x0, y0], [vx0, vy0]
    acc is a function defining the acceleration affecting
    the body to be simulated.
    returns x (orbital position) and v (orbital velocity)
    '''
    x, v = nt.leapfrog(x0, v0, t, acc)
    return x, v

def system_gravity(m1, m2, x):
    x = x.transpose()
    return (-vars.G*m1*m2/nt.norm(x)**3*x).transpose()

def system_acceleration(m, r, index, n):
    mask_index = np.arange(n) != index
    r_planet = r[index]
    r_between = r[index] - r[mask_index]
    acc = system_gravity(m[mask_index], m[index], r_between)/m[index]
    acc = np.sum(acc, axis = 0)
    return acc

def center_of_mass(m, r):
    total_mass = np.sum(m)
    cm = np.sum(m*r.transpose(), axis = 1)/total_mass
    return cm

def n_body_problem(xx, vv, cm, vcm, mass, time, n):
    v = np.zeros(vv[0].shape)
    dt = time[1] - time[0]
    for k in range(len(time)-1):
        x = np.copy(xx[k])
        for i in range(n):
            acc = lambda r: system_acceleration(mass, r, i, n)
            x[i] = xx[k,i] + vv[k,i]*dt + 0.5*acc(xx[k])*dt**2
        for i in range(n):
            acc = lambda r: system_acceleration(mass, r, i, n)
            v[i] = vv[k,i] + 0.5*(acc(xx[k])+ acc(x))*dt
        cm[k+1] = center_of_mass(mass, x)
        vcm[k] = (cm[k+1]-cm[k])/dt
        xx[k+1] = x
        vv[k+1] = v
    vcm[-1] = vcm[-2] # Approximate final value
    return xx, vv, cm, vcm

def n_body_setup(masses, time, steps, x0, v0, ref_frame = 'cm'):
    x0 = x0.transpose()
    v0 = v0.transpose()
    n = len(masses)
    cm  = np.zeros((steps, 2))
    vcm = np.zeros((steps, 2))
    x = np.zeros((len(time), n, 2))
    v = np.zeros((len(time), n, 2))
    x[0] = x0
    v[0] = v0
    xx, vv, cm, vcm = n_body_problem(x, v, cm, vcm, masses, time, n)
    xx = xx.transpose()
    vv = vv.transpose()

    if ref_frame == 'cm':
        for i in range(2):
            xx[i] = xx[i] - cm[:,i]
            vv[i] = vv[i] - vcm[:,i]
    return xx, vv, cm, vcm

import numpy as np
import numtools as nt
import variables as vars
from scipy.interpolate import interp1d

def gravity(m1, m2, x):
    return -vars.G*m1*m2/nt.norm(x)**3*x

def kepler3(m1, m2, a):
    return np.sqrt(4*np.pi**2/(vars.G*(m1+m2))*a**3)

def trajectory(masses, x, v, host, sat, target, sun, time, launched, tol):
    theta_target = np.arctan2(x[:, target, 1], x[:,target,0]) + np.pi
    theta_host = np.arctan2(x[:, host, 1], x[:,host,0]) + np.pi
    r_target = nt.norm(x[:, target], ax = 1)
    r_host = nt.norm(x[:, host], ax = 1)
    delta_v_peri = []
    launch_window = []
    t_cept = []
    semimajor = []
    for i in range(len(time)):
        r1 = r_host[i]
        t1 = theta_host[i]
        check = colinear(t1, theta_target, tol) # Check future values
        possibles = np.argwhere(check[i:] != 0) + i     # Values where planets align through sun.
        for possible in possibles:
            r2 = r_target[possible]
            a = (r1 + r2)/2
            p = kepler3(masses[sun], masses[sat], a)
            t_future = time[i] + p/2
            i_future = np.argmin(np.abs(time-t_future))
            if i_future <= len(time):
                if nt.norm(x[i_future, host] - x[possible, host], ax = 1) < tol: #i_future == possible
                    print('Found possible launch window')
                    transfer_peri = vis_viva(masses[sun], r1, a)*nt.unit_vector(v[i, host])
                    v_soi = transfer_peri - v[i, host]
                    if launched == True:
                        v_escape = 0
                    else:
                        v_escape = vis_viva(masses[host], vars.radius[0]*1000/vars.AU_tall, 1e20)
                    vfin = np.sqrt(v_escape**2 + nt.norm(v_soi)**2)
                    delta_v_peri.append(vfin)
                    launch_window.append(time[i])
                    t_cept.append(time[i_future])
                    semimajor.append(a)

    return delta_v_peri, launch_window, t_cept, semimajor

def trajectory_interp(masses, x, v, host, sat, target, sun, time, launched, tol):
    theta_target = np.arctan2(x[:, target, 1], x[:,target,0]) + np.pi
    theta_host = np.arctan2(x[:, host, 1], x[:,host,0]) + np.pi
    r_target = nt.norm(x[:, target], ax = 1)
    r_host = nt.norm(x[:, host], ax = 1)
    switch_target = np.argwhere(np.abs(np.diff(theta_target)) > np.pi)[0]
    switch_host = np.argwhere(np.abs(np.diff(theta_host)) > np.pi)[0]
    for i in switch_target:
        theta_target[i+1:] = theta_target[i+1:] + 2*np.pi
    for j in switch_host:
        theta_host[j+1:] = theta_target[j+1:] + 2*np.pi
    print(theta_target, theta_host)

    r_target = interp1d(theta_target, r_target)
    x = x.transpose()
    x_host = nt.interp_xin(time, x[:, host])
    x_target_t = nt.interp_xin(time, x[:, target])
    x_target_t = lambda t: np.array([x_target_t[0](t),  x_target_t[1](t)])
    x_host = lambda t: np.array([x_host[0](t), x_host[1](t)])

    x_theta = nt.interp_xin(theta_target, x[:, target])
    x_target_theta = lambda theta: np.array([x_theta[0](theta), x_theta[1](theta)])

    x = x.transpose()

    delta_v_peri = []
    launch_window = []
    t_cept = []
    semimajor = []
    for i in range(len(x[:, host])):
        r1 = r_host[i]
        match = colinear(theta_host[i], theta_target_t(t_future))
        possible_t = theta_target_t()
        if match:
            a = (r1 + r2)/2
            p = kepler3(masses[sun], masses[sat], a)
            t_future = time[i] + p/2

    return delta_v_peri, launch_window, t_cept, semimajor

def sphere_of_influence(a, m1, m2):
    return a*(m1/m2)**(2/5)

def colinear(theta1, theta2, tol):
    return np.abs(np.abs(theta1 - theta2) - np.pi) <= tol

def vis_viva(m_senter, r, a):
    return np.sqrt(vars.G*m_senter*(2/r - 1/a))

def grav_influence(m_star, m_planet, r_to_star, k = 10):
    return nt.norm(r_to_star, ax = 1)*np.sqrt(m_planet/(k*m_star))

def circularize(m_senter, m_sat, x, v, desired_r):
    '''
    m_senter is mass of center body in two body system
    m_sat is mass of body orbiting central body in an elliptical orbit
    x is the position of the orbiting body
    v is the velocity of the orbiting body
    desired_r is the radius of the desired circular orbit
    tol is the allowed deviation from the circular radius.
    '''
    v_desired = vis_viva(m_senter, desired_r, desired_r)
    delta_v = 0 #foob
    dv = np.zeros(x.shape)
    apoapsis = np.argmax(nt.norm(x, ax = 1))
    dv[apoapsis] = 0 #baz


def orbit(x0, v0, acc, t):
    '''
    Assumes x0, v0, are arrays containing [x0, y0], [vx0, vy0]
    acc is a function defining the acceleration affecting
    the body to be simulated.
    returns x (orbital position) and v (orbital velocity)
    '''
    x, v = nt.leapfrog(x0, v0, t, acc)
    return x, v

def patched_conic_orbits(time, mass, x0, v0, ref_frame = 0):
    '''
    time is time array
    ref_frame is the index of the planet/star of the desired reference
    frame, default is sun (index 0).
    Assumes mass is array of masses, with index 0 corresponding to sun
    x0 is array of initial positions, assumed to be of shape
    (1, number of planets, 2).
    v0 is array of initial velocities, same shape.
    '''
    steps = len(time)
    x = np.zeros((steps, x0.shape[0], x0.shape[1]))
    v = np.zeros(x.shape)
    for i in range(1, len(mass)):
        acc = lambda r, t: gravity(mass[0], mass[i], r)/mass[i]
        a, b = orbit(x0[i], v0[i], acc, time)
        x[:, i], v[:, i] = a.transpose(), b.transpose()
        x[:, i] = x[:,i] - x[:,ref_frame]
        v[:, i] = v[:,i] - v[:,ref_frame]
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
    for k in range(len(time)-1):
        dt = time[k+1] - time[k]
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

def n_body_sat(xp, mass, time, dv, sx0, sv0, sm, opt_vel = None, opt_orb = None, t_opt = 0, opt = False, numboosts = 1000):
    def acc(r_sat, r):
        r_between = r_sat - r
        rr1 = nt.norm(r_between, ax = 1)
        acc = 0
        for mm, rr, rb in zip(mass, rr1, r_between):
            acc1 = -vars.G*mm/rr**3*rb
            acc += acc1
        return acc

    x_sat = np.zeros((len(time), 2))
    v_sat = np.zeros((len(time), 2))
    x_sat[0] = sx0
    v_sat[0] = sv0

    nextboost = np.argmin(np.abs(time-t_opt))

    for k in range(len(time) -1):
        dt = time[k+1] - time[k]
        acc1 = acc(x_sat[k], xp[k])
        x_sat[k+1] = x_sat[k] + v_sat[k]*dt + 0.5*acc1*dt**2
        acc2 = acc(x_sat[k+1], xp[k+1])
        if opt and time[k] > t_opt and k == nextboost:
            v_diff = opt_vel[k + 1] - v_sat[k]
            if nt.norm(v_diff) < 0.1:
                dv[k] = dv[k] + v_diff
            v_sat[k+1] = v_sat[k] + 0.5*(acc1 + acc2)*dt + dv[k]
            nextboost += int(len(time)/numboosts)
        else:
            v_sat[k+1] = v_sat[k] + 0.5*(acc1 + acc2)*dt + dv[k]
    return x_sat, v_sat, dv

def n_body_setup(masses, time, steps, x0, v0, ref_frame = 'cm'):
    #x0 = x0.transpose()
    #v0 = v0.transpose()
    n = len(masses)
    cm  = np.zeros((steps, 2))
    vcm = np.zeros((steps, 2))
    x = np.zeros((len(time), n, 2))
    v = np.zeros((len(time), n, 2))
    x[0] = x0
    v[0] = v0
    cm[0] = center_of_mass(masses, x0)
    xx, vv, cm, vcm = n_body_problem(x, v, cm, vcm, masses, time, n)
    xx = xx.transpose()
    vv = vv.transpose()
    if ref_frame == 'cm':
        for i in range(2):
            xx[i] = xx[i] - cm[:,i]
            vv[i] = vv[i] - vcm[:,i]
    elif ref_frame == 'sol':
        for i in range(2):
            xx[i] = xx[i] - xx[i,0]
            vv[i] = vv[i] - vv[i,0]
    return xx, vv, cm, vcm

import numpy as np
import numtools as nt
import variables as vars

def gravity(m1, m2, x):
    return -vars.G*m1*m2/nt.norm(x)**3*x

def kepler3(m1, m2, a):
    return np.sqrt(4*np.pi**2/(vars.G*(m1+m2))*a**3)

def trajectory(masses, x, v, steps, host, sat, target, sun, time, launched, tol, i_tol = 100):
    theta_target = np.arctan2(x[:, target, 1], x[:,target,0]) + np.pi
    theta_host = np.arctan2(x[:, host, 1], x[:,host,0]) + np.pi

    r_target = nt.norm(x[:, target], ax = 1)
    r_host = nt.norm(x[:, host], ax = 1)

    delta_v_peri = []
    launch_window = []

    dt = time[1] - time[0]

    for i in range(len(x[:, host])):
        r1 = r_host[i]
        t1 = theta_host[i]
        check = colinear(t1, theta_target, tol) # Check future values
        possibles = np.argwhere(check[i:] != 0) + i     # Values where planets align through sun.
        for possible in possibles:
            r2 = r_target[possible]
            a = (r1 + r2)/2
            p = kepler3(masses[sun], masses[sat], a)
            i_future =  int(p/(2*dt))     #index of intercept = i + i_future
            if abs(i + i_future - possible) < i_tol:
                try:
                    if nt.norm(x[i + i_future, target] - x[possible, target], ax = 1) < tol:
                        transfer_peri = vis_viva(masses[sun], r1, a)*nt.unit_vector(v[i, host])
                        v_soi = transfer_peri - v[i, host]
                        if launched == True:
                            v_escape = 0
                        else:
                            v_escape = vis_viva(masses[host], vars.radius[0]*1000/vars.AU_tall, 1e20)
                        vfin = np.sqrt(v_escape**2 + nt.norm(v_soi)**2)
                        delta_v_peri.append(vfin)
                        launch_window.append(time[i])
                except IndexError:
                    break
    return delta_v_peri, launch_window

def sphere_of_influence(a, m1, m2):
    return a*(m1/m2)**(2/5)

def colinear(theta1, theta2, tol):
    return np.abs(np.abs(theta1 - theta2) - np.pi) < tol

def vis_viva(m_senter, r, a):
    return np.sqrt(vars.G*m_senter*(2/r - 1/a))

def simple_potential(radii, masses, body_index):
    ep = 0
    for radius, mass in radii, masses:
        ep -= vars.G*masses[body_index]*mass*(1/radius[1] - 1/radius[0])
    return ep

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

def n_body_sat(xp, vp, mass, time, host, dv, launched, sx0, sv0, sm):
    def acc(r_sat, r):
        r_between = r_sat - r
        acc = system_gravity(mass, sm, r_between)/sm
        acc = np.sum(acc, axis = 0)
        return acc

    x_sat = np.zeros((len(time), 2))
    v_sat = np.zeros((len(time), 2))
    x_sat[0] = sx0
    v_sat[0] = sv0

    dt = time[1] - time[0]
    for k in range(len(time) -1):

        if launched == False:
            x_sat[k] = xp[k, host] + vars.radius[0]*1000/vars.AU_tall*nt.unit_vector(vp[k,host])
            v_sat[k] = vp[k, host]

        if launched == False and dv[k] != 0:
            x_sat[k] = xp[k, host] + vars.radius[0]*1000/vars.AU_tall*nt.unit_vector(vp[k,host])
            v_sat[k] = vp[k, host] + dv[k]*nt.unit_vector(vp[k, host])
            launched = True

        if launched == True:
            acc1 = acc(x_sat[k], xp[k])
            x_sat[k+1] = x_sat[k] + v_sat[k]*dt + 0.5*acc1*dt**2
            acc2 = acc(x_sat[k+1], xp[k+1])
            v_sat[k+1] = v_sat[k] + 0.5*(acc1 + acc2)*dt

        if launched == False:
            x_sat[k] = xp[k, host] + vars.radius[0]*1000/vars.AU_tall*nt.unit_vector(vp[k,host])
            v_sat[k] = vp[k, host]

        if launched == False and dv[k] != 0:
            x_sat[k] = xp[k, host] + vars.radius[0]*1000/vars.AU_tall*nt.unit_vector(vp[k,host])
            v_sat[k] = vp[k, host] + dv[k]*nt.unit_vector(vp[k, host])
            launched = True

    return x_sat, v_sat



def n_body_custom(mass, t, x, v, host, dv, launched, sx0, sv0, sm):
    xx, vv = n_body_sat(x, v, mass, t, host, dv, launched, sx0, sv0, sm)
    xx = xx.transpose()
    vv = vv.transpose()
    return xx, vv


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
    return xx, vv, cm, vcm

def orbital_maneuver(system_x0, system_v0, sat_x0, sat_v0, target_index):
    pass

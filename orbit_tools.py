import numpy as np
import numtools as nt
import variables as vars

def gravity(m1, m2, x):
    return -vars.G*m1*m2/nt.norm(x)**3*x

def kepler3(m1, m2, a):
    return np.sqrt(4*np.pi**2/(vars.G*(m1+m2))*a**3)

def trajectory(system_masses, system_x, steps, host, sat, target, sun, dt):

    theta_target = np.arctan2(system_x[:, target, 1], system_x[:,target,0]) + np.pi
    theta_host = np.arctan2(system_x[:, host, 0], system_x[:,host,1]) + np.pi
    r_target = nt.norm(system_x[:, target], ax = 1)
    r_host = nt.norm(system_x[:, host], ax = 1)
    tol = 0.000
    import matplotlib.pyplot as plt

    for i in range(len(system_x[:, host])):
        r1 = r_host[i]
        t1 = theta_host[i]
        check = np.where(abs(np.abs(t1 - theta_target) - np.pi) < tol, theta_target, 0)
        possibles = np.argwhere(check != 0)
        for possible in possibles:
            r2 = r_target[possible]
            a = (r1 + r2)/2
            p = kepler3(system_masses[sun], system_masses[sat], a)
            #print(int(p/(2*dt)))
            try:
                if nt.norm(system_x[i + int(p/(2*dt)), target] - system_x[possible, target], ax = 1) < tol:
                    print(i*dt)
            except IndexError:
                break
                print('Index out of range')

def simple_potential(radii, masses, body_index):
    ep = 0
    for radius, mass in radii, masses:
        ep -= vars.G*masses[body_index]*mass*(1/radius[1] - 1/radius[0])
    return ep

    #approx = vars.G*(system_masses[body_index]*vars.m_star)/(np.sqrt(x0**2+y0**2)) - vars.G*(system_masses[body_index]*vars.m_star)/(np.sqrt(x1**2+y1**2))
    #print(approx)
    # AU = m/AUtall
    # år = s/faktor
    #solm = kg/solmas
    # E = kgm**2/s = (solm*solmass)*(AU*AUtall)**2/(år*faktor)**2


    #e_p = -(np.trapz(fx, x) + np.trapz(fy, y))



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

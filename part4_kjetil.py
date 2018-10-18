import numpy as np
import variables as vars
from numpy.linalg import inv as numpy_inv
from scipy.interpolate import interp1d as numpy_interpolate

def vel_rel_star(Dlam, lam0):
    vr = Dlam/lam0*vars.c #nm/nm*m/s -> [m/s]
    return vr[:,0] + vr[:,1]

def shift_ref(p, d, convert = 'from'):
    '''convert 'to' or 'from' the given coordinate system'''
    if convert == 'from':
        M = np.array([np.cos(p), np.sin(p)]).transpose()
        M_inv = numpy_inv(M)
        return M_inv * d
    elif convert == 'to':
        M = np.array([np.cos(p), np.sin(p)]).transpose()
        return  M * d
    else:
        print('convert has to be either to, or from. Default is from')
        pass

def velocity_from_stars(lam_measured, lam0 = 656.3): #phi [rad], lam [nm]
    '''calculates velocity in cartesian coordinates given angles and delta lambdas'''
    ref_stars = np.array(vars.ref_stars) #phi in DEGREES, lam in nanometers
    phi_skew = np.array(ref_stars[:,0])*np.pi/180 #phi in RADIANS
    lam_skew = np.array(ref_stars[:,1])
    lam_delta = lam_skew - lam_measured #difference between reference shift and measured shift

    lam = shift_ref(phi_skew, lam_delta, 'from')
    vel_cart = vel_rel_star(lam, lam0)# [m/s]
    #print('velocity of satelite with respect to sun in cartesian coordinates', vel_cart)
    return vel_cart/vars.AU_tall*vars.year

#TRILATERATION
#time of measurement: t0
#list of meadured distances: [p0, p1... pn, star]

def position_from_objects(current_time, distances, xx):
    #print('DIST', distances)
    d = np.zeros(len(distances))
    d[0] = distances[-1]
    d[1:] = distances[:-1]
    #print('d', d)
    #xx = np.load('saved/saved_orbits/launch_resolution/pos.npy')
    p = xx[:,:,current_time].transpose()
    #print('p', p)
    #d = np.random.random(9)     #distances
    #p = np.random.random([9,2]) #positions

    x = np.zeros([])
    y = np.zeros([])
    count = 0
    ind = np.linspace(0,len(d)-1,len(d), dtype = 'int')
    for i in ind:
        for j in ind[ind != i]:
            for k in ind[(ind != i) * (ind != j)]:
                if i != 0:
                    #print(i,j,k)

                    count += 1
                    a2 = p[j,0] - p[i,0]    #a corresponds to x positions of planets
                    a3 = p[k,0] - p[i,0]
                    b2 = p[j,1] - p[i,1]    #b corresponds to y positions of planets
                    b3 = p[k,1] - p[i,1]
                    #c is a constant depandant on x and y positions in addition to distances
                    c2 = d[i]**2 - d[j]**2 - p[i,0]**2 - p[i,1]**2 + p[j,0]**2 + p[j,1]**2
                    c3 = d[i]**2 - d[k]**2 - p[i,0]**2 - p[i,1]**2 + p[k,0]**2 + p[k,1]**2

                    y = np.append(y, 1/2*(a2*c3 - a3*c2) / (a2*b3 - a3*b2))
                    x = np.append(x, (c3 - 2*y[i]*b3) / (2*a3))
                else:
                    pass
    #print(x)
    x_avg = np.average(x)
    y_avg = np.average(y)
    pos = np.array([x_avg, y_avg])
    return pos

if __name__ == '__main__':
    #lam_measured = np.array([0.00128*2, -0.005186*2]) #measured delta_lambda
    #velocity_from_stars(lam_measured)
    x, y = position_from_objects(time_of_measurement, list_of_measured_distances)

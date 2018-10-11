import numpy as np
import variables as vars
from numpy.linalg import inv as numpy_inv

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
    lam_delta = lam_measured - lam_skew #

    lam = shift_ref(phi_skew, lam_delta, 'from')
    vel_cart = vel_rel_star(lam, lam0)# [m/s]
    print('velocity of satelite with respect to sun in cartesian coordinates', vel_cart)

#TRILATERATION
#time of measurement: t0
#list of meadured distances: [p0, p1... pn, star]

def position_from_objects(current_time, distances):
    #distances =  SOLARSYSTEM.analyse_distances
    #positions = xx[:,:,current_time]
    d = np.random.random(9)     #distances
    p = np.random.random([9,2]) #positions
    x = np.zeros(len(d)-2)
    y = np.zeros(len(d)-2)

    for i in range(len(d)-2):
        #defining constants to make the final expression readable
        a2 = p[i+1,0] - p[i,0]    #a corresponds to x positions of planets
        a3 = p[i+2,0] - p[i,0]
        b2 = p[i+1,1] - p[i,1]    #b corresponds to y positions of planets
        b3 = p[i+2,1] - p[i,1]
        #c is a constant depandant on x and y positions in addition to distances
        c2 = d[i]**2 - d[i+1]**2 - p[i,0]**2 - p[i,1]**2 + p[i+1,0]**2 + p[i+1,1]**2
        c3 = d[i]**2 - d[i+2]**2 - p[i,0]**2 - p[i,1]**2 + p[i+2,0]**2 + p[i+2,1]**2

        y[i] = 1/2*(a2*c3 - a3*c2) / (a2*b3 - a3*b2)
        x[i] = (c3 - 2*y[i]*b3) / (2*a3)
    print(x)
    print(y)
    x_avg = np.average(x)
    y_avg = np.average(y)
    print(x_avg, y_avg)
    return x_avg, y_avg

if __name__ == '__main__':
    #lam_measured = np.array([0.00128*2, -0.005186*2]) #measured delta_lambda
    #velocity_from_stars(lam_measured)
    x, y = position_from_objects(time_of_measurement, list_of_measured_distances)

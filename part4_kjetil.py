import numpy as np
import variables as vars
#dopplershift analysis

def vel_rel_star(Dlam1, Dlam2, lam0 = 656.3):
    vr1 = Dlam1/lam0*vars.c #nm/nm*m/s -> [m/s]
    vr2 = Dlam2/lam0*vars.c
    return vr1, vr2
#vr = Dlam/lam0*c
#def vel_rel_star():
star1 = np.array(vars.ref_stars[0])
star1[0] = star1[0]*2*np.pi/360
star2 = np.array(vars.ref_stars[1])
star2[0] = star2[0]*2*np.pi/360
lam0 = 656.3 #nm
print(star1) #phi in radians, lamb in nanometers
print(star2) #phi in radians, lamb in nanometers

vel1_p2star, vel2_p2star = vel_rel_star(star1[1], star2[1]) #m/s
vel1_s2star, vel2_s2star = vel_rel_star(0, 0)
vel1 = vel1_p2star - vel1_s2star
vel2 = vel2_p2star - vel2_s2star
print(vel1,vel2) #wrong referance system
#def shift_ref(x, y, th1, th2):

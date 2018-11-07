import numpy as np
import matplotlib.pyplot as plt
import variables as vars
#import part6_kjetil_density as part6_kd
import part6_kjetil as part6
import numtools as nt

def new_coordinates(old, t): #polar coordinates [r, theta]
    return old[1] + 2*np.pi/period*t

def c2p_pos(cart): #x, y
    polar = np.zeros(cart.shape)
    polar[0] = nt.norm(cart) #numtools
    polar[1] = np.arctan2(cart[1], cart[0])
    return polar

def p2c_pos(polar): #r, theta #USEFUL
    cart = np.zeros(polar.shape)
    cart[0] = polar[0]*np.cos(polar[1])
    cart[1] = polar[0]*np.sin(polar[1])
    return cart

def c2p_vel(pos_cart, vel_cart): #x, y, v_x, v_y
    polar = np.zeros(pos_cart.shape)
    polar[0] = (pos_cart[0]*vel_cart[0] + pos_cart[1]*vel_cart[1])/nt.norm(pos_cart) #numtools
    polar[1] = (pos_cart[0]*vel_cart[1] - vel_cart[0]*pos_cart[1])/(pos_cart[0]**2 + pos_cart[1]**2)
    return polar

def p2c_vel(pos_polar, vel_polar): #r, theta ??????????????????++
    cart = np.zeros(pos_polar.shape)
    cart[0] = vel_polar[0]*np.cos(pos_polar[1]) - pos_polar[0]*np.sin(pos_polar[1])*vel_polar[1]
    cart[1] = vel_polar[0]*np.sin(pos_polar[1]) + pos_polar[0]*np.cos(pos_polar[1])*vel_polar[1]
    return cart



def landing(pos, vel, boost_time, boost): #Where pos and vel is given relative to the planet, not relative to the sun
    #pos = c2t_pos(pos_c) #position relative to planet in tangential coordinates
    #vel = c2t_vel(pos_c, vel_c) #velocity relative to planet in tangential coordinates
    period = vars.period[1]*24*60*60
    radius = vars.radius_normal_unit[1]
    A = vars.area_lander  #cross sectional area
    M_planet = vars.m_normal_unit[1]
    m = vars.mass_lander
    G = vars.G_SI
    #r_vec = np.array([1,0])
    polar_vec = np.array([0,1])
    dt = 0.1
    t = 0
    count = 0
    boosted = False
    #theta is not really theta, but the tangent given in meters
    while nt.norm(pos) > radius:
    #while count < 80000:
        wel = p2c_vel(c2p_pos(pos), 2*np.pi/period * polar_vec)
        vel_drag = vel - wel
        rho = part6.density(nt.norm(pos))
        acc = -1/2*rho*A*nt.norm(vel_drag)*vel_drag/m - (M_planet*G/(nt.norm(pos)**2))*nt.unit_vector(pos)
        vel = vel + acc*dt
        pos = pos + vel*dt
        t += dt
        count += 1
        if t >= boost_time:
            if not boosted:
                plt.scatter(pos[0], pos[1], c = 'b')
                boosted = True
                angle1 = c2p_pos(pos)[1]
                vel = vel*boost
        if count % int(50) == 0:
            plt.scatter(pos[0], pos[1], 1, 'r')
    angle2 = c2p_pos(pos)[1]
    plt.scatter(pos[0], pos[1], c = 'b')
    print('DRAG: You reached the surface with a velocity of %.3f m/s after %.2f hours' %(nt.norm(vel_drag), (t-boost_time)/60/60))

    vel_drag_radial = nt.norm(vel_drag*nt.unit_vector(-pos))
    print(vel_drag_radial)
    vel_drag_tangential = nt.norm(vel_drag*nt.rotate(nt.unit_vector(-pos), np.pi/2))
    print(vel_drag_tangential)

    print('DRAG: Radial velocity was %.3f m/s and tangential velocity was %.3f m/s' %(vel_drag_radial, vel_drag_tangential))
    print('Angle', angle2-angle1)
    plt.axis('equal')
    pi_vec = np.linspace(0, 2*np.pi, 1000)
    for theta in pi_vec:
        circle = p2c_pos(np.array([radius, theta]))
        plt.scatter(circle[0], circle[1], 0.1, 'k')
    plt.show()
    return angle1, angle2-angle1, boost_time - t, pos, vel_drag

if __name__ == '__main__':
    position = np.array([vars.radius_normal_unit[1] + 400000, 0])
    velocity = np.array([0, 2750])
    boost = 0.8
    boost_time = 1000
    alpha, beta, duration, pos, vel_drag = landing(position, velocity, boost_time, boost)
    

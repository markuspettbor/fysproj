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



def landing(pos, vel, boost_angle, boost, plot = False): #Where pos and vel is given relative to the planet, not relative to the sun
    plt.figure()
    #pos = c2t_pos(pos_c) #position relative to planet in tangential coordinates
    #vel = c2t_vel(pos_c, vel_c) #velocity relative to planet in tangential coordinates
    period = vars.period[1]*24*60*60
    radius = vars.radius_normal_unit[1]
    A = vars.area_lander  #cross sectional area
    A_parachute = 25
    M_planet = vars.m_normal_unit[1]
    m = vars.mass_lander
    G = vars.G_SI
    #r_vec = np.array([1,0])
    polar_vec = np.array([0,1])
    dt = 0.1
    t = 0
    count = 0
    boosted = False
    parachuted = False
    angle_less = False
    #print('POS IN LANDING', pos)
    #print('VEL IN LANDING', vel)
    parachute_time = 0
    #theta is not really theta, but the tangent given in meters
    while nt.norm(pos) > radius:
    #while count < 140000:
        wel = p2c_vel(c2p_pos(pos), 2*np.pi/period * polar_vec)
        vel_drag = vel - wel
        rho = part6.density(nt.norm(pos))
        acc = -1/2*rho*A*nt.norm(vel_drag)*vel_drag/m - (M_planet*G/(nt.norm(pos)**2))*nt.unit_vector(pos)
        vel = vel + acc*dt
        pos = pos + vel*dt
        t += dt
        count += 1
        #print(boost_angle , 'POSS')
        if np.arctan2(pos[1],pos[0]) >= boost_angle and angle_less:
            if not boosted:
                if plot:
                    plt.scatter(pos[0], pos[1], c = 'b')
                print('boosted at angle', np.arctan2(pos[1],pos[0]))
                print('boost angle is  ', boost_angle)
                boosted = True
                boost_time = t
                print('boost_time', boost_time)
                angle1 = np.arctan2(pos[1],pos[0])
                vel = vel*boost
        elif np.arctan2(pos[1],pos[0]) < boost_angle:
            angle_less = True
        if nt.norm(vel_drag) < 100 and not parachuted:
            parachute_time = t
            parachuted = True
            print('paratime', parachute_time)
            A = A_parachute
            acc = -1/2*rho*A*nt.norm(vel_drag)*vel_drag/m - (M_planet*G/(nt.norm(pos)**2))*nt.unit_vector(pos)
            F = m*acc
            print('Force', nt.norm(F), 'N')
        if count % int(50) == 0 and plot:
            plt.scatter(pos[0], pos[1], 1, 'r')
    angle2 = np.arctan2(pos[1],pos[0])
    if plot:
        print('time', t)
        plt.scatter(pos[0], pos[1], c = 'b')
        print('DRAG: You reached the surface with a velocity of %.3f m/s after %.2f seconds' %(nt.norm(vel_drag), (t-boost_time)))

        vel_drag_radial = nt.norm(vel_drag*nt.unit_vector(-pos))
        #print(vel_drag_radial)
        vel_drag_tangential = nt.norm(vel_drag*nt.rotate(nt.unit_vector(-pos), np.pi/2))
        #print(vel_drag_tangential)

        print('DRAG: Radial velocity was %.3f m/s and tangential velocity was %.3f m/s' %(vel_drag_radial, vel_drag_tangential))
        print('Angle', angle2-angle1)
        plt.axis('equal')
        pi_vec = np.linspace(0, 2*np.pi, 1000)
        for theta in pi_vec:
            circle = p2c_pos(np.array([radius, theta]))
            circle2 = circle * 1.27
            plt.scatter(circle[0], circle[1], 0.1, 'k')
            plt.scatter(circle2[0], circle2[1], 0.1, 'k')
        #plt.show()
    return angle2, angle2-angle1, parachute_time, boost_time, pos, vel_drag

def optimise_landing(pos_o, vel_o, angle_landing, boost, plotting = False):
    period = vars.period[1]*24*60*60
    radius = vars.radius_normal_unit[1]
    pos_o = pos_o*vars.AU_tall
    vel_o = vel_o*vars.AU_tall/vars.year
    #angle_initial, alpha_initial, time_initial = landing(position_vec[start_index], velocity_vec[start_index], angle_landing, boost, plot = True)[0:3]
    accuracy = 3
    angle_last = angle_landing
    beta = 0
    for i in range(accuracy):
        angle_release = angle_last - beta
        angle, alpha, time_para, time_boost = landing(pos_o, vel_o, angle_release, boost, plot = plotting)[0:4]
        angle_last = angle_release
        beta = angle - angle_landing
        print('error', beta)
        #w = 2*pi/period
    return time_para, time_boost

if __name__ == '__main__':
    time, position, velocity, position_planet, angle, index = np.load('saved/saved_orbits/data_to_lander.npy')
    time = time*vars.year
    position = position*vars.AU_tall
    velocity = velocity*vars.AU_tall/vars.year
    position_planet = position_planet*vars.AU_tall
    velocity_planet = np.gradient(position_planet, axis = 0)/(time[1]-time[0])
    #print((position_planet[1]- position_planet[0])/(time[1]-time[0]))
    #print(velocity_planet)
    #print(velocity)
    velocity = (velocity - velocity_planet)
    #print(velocity)

    #position = np.array([vars.radius_normal_unit[1] + 400velocity[-1,:]000, 0])
    #velocity = np.array([0, 2750])
    boost = 0.8
    boost_time = 1000
    '''
    plt.plot(position[:,0], position[:,1])
    plt.scatter(position[325,0], position[325,1], c = 'r')
    plt.scatter(position[350,0], position[350,1], c = 'g')
    plt.scatter(position[375,0], position[375,1], c = 'b')
    plt.scatter(position[400,0], position[400,1], c = 'r')
    plt.scatter(position[425,0], position[425,1], c = 'g')
    plt.scatter(position[450,0], position[450,1], c = 'b')
    plt.scatter(position[475,0], position[475,1], c = 'r')
    plt.axis('equal')
    plt.show()
    '''
    index = -2

    print('printy', position[index,:])
    #alpha, beta, duration, pos, vel_drag = landing(position[index], velocity[index], angle, boost, plot = True)
    release_angle, time_para = optimise_landing(position[index,:], velocity[index,:], angle, boost, plotting = True)
    plt.show()
    plt.plot(position[:,0], position[:,1])
    plt.show()
    print('Release angle:', release_angle)
    print('Time parachute:', time_para)

import numpy as np
import matplotlib.pyplot as plt
import variables as vars
import part6
import numtools as nt

# This code is ours

def new_angle(angle_old, t_measured, t_landing): #polar coordinates [r, theta]
    t = t_landing - t_measured
    period = vars.period[1]*24*60*60 #sec
    angle_new = angle_old + 2*np.pi/period*t
    return angle_new

def eject_time(t_measured, wp, ws, beta0, beta, alpha):
    '''t_measured is time angle is measured
    wp and ws is angular velocity of planet and satellite
    beta0 is angle of satellite at t = 0
    beta is angle it takes to land (could vary a little)
    alpha is angle of landingsite at measurement time'''
    return (beta0 + beta - alpha + wp*t_measured)/(wp - ws)

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
    #plt.figure()
    period = vars.period[1]*24*60*60
    radius = vars.radius_normal_unit[1]
    print('RADIUS[1]', radius)
    A = vars.area_lander  #cross sectional area
    A_parachute = 25
    M_planet = vars.m_normal_unit[1]
    m = vars.mass_lander
    G = vars.G_SI
    polar_vec = np.array([0,1])
    dt = 0.1
    t = 0
    count = 0
    boost_time = 1000
    parachute_time = 1001
    angle1 = 0
    boost_velocity = 0
    time_force = np.linspace(0, 30000, 300001)
    force = np.zeros(300001)
    print('forceShapes', time_force.shape, force.shape)
    boosted = False
    parachuted = False
    angle_less = False
    #theta is not really theta, but the tangent given in meters
    print('NORM', nt.norm(pos))

    while abs(boost_angle) > np.pi:
        if boost_angle > 0:
            boost_angle = boost_angle - 2*np.pi
        elif boost_angle < 0:
            boost_angle = boost_angle + 2*np.pi

    while nt.norm(pos) > radius:
    #while count < 240000:
        wel = p2c_vel(c2p_pos(pos), 2*np.pi/period * polar_vec)
        vel_drag = vel - wel
        rho = part6.density(nt.norm(pos))
        acc = -1/2*rho*A*nt.norm(vel_drag)*vel_drag/m - (M_planet*G/(nt.norm(pos)**2))*nt.unit_vector(pos)
        vel = vel + acc*dt
        pos = pos + vel*dt
        t += dt
        force[count] = nt.norm(acc*m)
        count += 1
        if np.arctan2(pos[1],pos[0]) >= boost_angle and angle_less:
            #a check that should trigger when the desired angle is reached, and not before this
            #without the 'and angle_less' it would not work because angles are stupid 0->2pi->0
            if not boosted:
                if plot:
                    plt.scatter(pos[0]/1000, pos[1]/1000, c = 'b')
                print('boosted at angle', np.arctan2(pos[1],pos[0]))
                print('boost angle is  ', boost_angle)
                boosted = True
                boost_time = t
                print('boost_time', boost_time)
                angle1 = np.arctan2(pos[1],pos[0])
                boost_velocity = -vel*(1-boost)
                vel = vel*boost

        elif np.arctan2(pos[1],pos[0]) <= boost_angle:
            angle_less = True
        if nt.norm(vel_drag) < 30 and not parachuted:
            #parachute will be launched when the drag velocity is below 30m/s
            #could have included a check to see if the new force will be above
            #some threshold, though this should not happen at 30m/s
            parachute_time = t
            parachuted = True
            print('paratime', parachute_time)
            A = A_parachute
            acc = -1/2*rho*A*nt.norm(vel_drag)*vel_drag/m - (M_planet*G/(nt.norm(pos)**2))*nt.unit_vector(pos)
            F = m*acc
            force[count] = nt.norm(acc*m)
            print('Force', nt.norm(F), 'N')
        if count % int(50) == 0 and plot:
            plt.scatter(pos[0]/1000, pos[1]/1000, 1, 'r')
    angle2 = np.arctan2(pos[1],pos[0])
    if plot:
        print('time', t)
        plt.scatter(pos[0]/1000, pos[1]/1000, c = 'b')
        print('DRAG: You reached the surface with a velocity of %.3f m/s after %.2f seconds' %(nt.norm(vel_drag), (t-boost_time)))
        vel_drag_radial = nt.norm(vel_drag*nt.unit_vector(-pos))
        vel_drag_tangential = nt.norm(vel_drag*nt.rotate(nt.unit_vector(-pos), np.pi/2))
        print('DRAG: Radial velocity was %.3f m/s and tangential velocity was %.3f m/s' %(vel_drag_radial, vel_drag_tangential))
        print('Angle', angle2-angle1)
        plt.axis('equal')
        plt.xlabel('X-position [km]')
        plt.ylabel('Y-position [km]')
        pi_vec = np.linspace(0, 2*np.pi, 1000)
        for theta in pi_vec:
            circle = p2c_pos(np.array([radius, theta]))
            #circle2 = circle * 1.27
            plt.scatter(circle[0]/1000, circle[1]/1000, 0.1, 'k')
            #plt.scatter(circle2[0], circle2[1], 0.1, 'k')
        plt.show()
        #plt.plot(time_force, force, '-k')
        #plt.xlabel('Time [s]')
        #plt.ylabel('Force [N]')
        #plt.show()
    print('BOOST TIME ', boost_time)
    return angle2, angle2-angle1, parachute_time, boost_time, boost_velocity, t

def optimise_landing(position, velocity, angle_landing, boost, plotting = False):
    period = vars.period[1]*24*60*60
    radius = vars.radius_normal_unit[1]
    accuracy = 4
    angle_last = angle_landing
    beta = 0
    for i in range(accuracy):
        angle_release = angle_last + beta
        angle, alpha, time_para, time_boost, boost_velocity, time_landed = landing(position, velocity , angle_release, boost, plot = plotting)
        angle_last = angle_release
        beta = new_angle(angle_landing, 6137, time_landed) - angle
        while abs(beta) > np.pi:
            if beta > 0:
                beta = beta - 2*np.pi
            elif beta < 0:
                beta = beta + 2*np.pi
        print('error', beta)
        print('ANGLE LANDED AFTER ROTATION', new_angle(angle, 6137, time_landed) - angle_landing)
    return time_para, time_boost, boost_velocity, time_landed, angle_last - beta

if __name__ == '__main__':
    position, velocity, angle, boost = np.load('saved/saved_orbits/data_l.npy')
    print(position)
    print(velocity)
    print(angle)
    print(boost)
    release_angle, time_para, boost_velocity, time_landed = optimise_landing(position, velocity, angle, boost, plotting = True)

    print('Release angle:', release_angle)
    print('Time parachute:', time_para)
    plt.show()

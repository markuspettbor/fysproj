import numpy as np
import matplotlib.pyplot as plt
import variables as vars
import numtools as nt

def Part9_Ex3_P1():
    dtau = 32.82
    ts = np.array([0, 33.36, 66.71, 1137, 1237, 1517])
    dts_first = ts[1]- ts[0]
    dts_last = (ts[-1] - ts[-2])
    vs = 0.202
    G = vars.G_SI
    c = vars.c
    sm = vars.solmasse
    AU = vars.AU_tall
    r = AU
    mbh = 1.24166e7
    mbh_meter = mbh*sm*G/c**2
    M_r = mbh_meter/r
    E_m = (1 - 2*mbh_meter/r)**(1/2)/(1-vs**2)**(1/2)
    print('black hole mass', mbh_meter)
    print('M/r', M_r)
    print('E/m', E_m)

    r_first = 2/(1-E_m*dtau/dts_first)
    r_last = 2/(1-E_m*dtau/dts_last)

    print('AU', AU )
    print('first', r_first)
    print('last', r_last)

    r_first_AU = 2*mbh_meter/(1-E_m*dtau/dts_first)
    r_last_AU = 2*mbh_meter/(1-E_m*dtau/dts_last)

    print('AU', 1 )
    print('first', r_first_AU/AU)
    print('last', r_last_AU/AU)


def Part9_Ex3_P2():
    signal_C1, time_C1 = np.loadtxt('part2C_5_frame1.txt')
    signal_C2, time_C2 = np.loadtxt('part2C_5_frame2.txt')
    signal_E1, time_E1 = np.loadtxt('part2E_2_frame1.txt')
    signal_E2, time_E2 = np.loadtxt('part2E_2_frame2.txt')

    plt.subplot(2,1,1)
    plt.title('Planet Frame')
    plt.plot(signal_C1[:-1], np.diff(time_C1), '.k', linewidth = 1)
    plt.plot(signal_E1[:-1], np.diff(time_E1), '+k', linewidth = 1)
    plt.xlabel('Signal')
    plt.ylabel('Time Difference [s]')
    plt.legend(['Event happens', 'Event is seen'])

    plt.subplot(2,1,2)
    plt.title('Spaceship Frame')
    plt.plot(signal_C2[:-1], np.diff(time_C2), '.k', linewidth = 1)
    plt.plot(signal_E2[:-1], np.diff(time_E2), '+k', linewidth = 1)
    plt.xlabel('Signal')
    plt.ylabel('Time Difference [s]')
    plt.legend(['Event happens', 'Event is seen'])

    plt.tight_layout()
    plt.show()

def Part9_Ex5():
    c = 299792458 #METER/SEKUND
    G = 6.67e-11    #konstanti
    mass = 3.459198972171e23 #KILOGRAM
    planet_radius = 2678881.4788 #METER
    t1 = 14.6189759
    ts1_t1 = 14.5934565
    ts2_t1 = 14.5898827
    t2 = 10247.9021048 #sek
    ts1_t2 = 10247.8821483
    ts2_t2 = 10247.8863525
    pos1_t1 = np.array([5396873, -4194131])    #METER
    pos2_t1 = np.array([6770895, -933788])     #METER
    pos1_t2 = np.array([-6586081, 1827708])    #METER
    pos2_t2 = np.array([-6617567, -1710199])     #METER
    radius1 = nt.norm(pos1_t1)
    radius2 = nt.norm(pos2_t1)
    radius = radius1/2 + radius2/2
    vel = np.sqrt(G*mass/radius) #METER PER SEC
    print('r1', radius1)
    print('r2', radius2)
    print('vel', vel)

    beta = np.array([np.arctan2(pos1_t1[1], pos1_t1[0]) + 2*np.pi, np.arctan2(pos2_t1[1], pos2_t1[0]) + 2*np.pi, \
                     np.arctan2(pos1_t2[1], pos1_t2[0]), np.arctan2(pos2_t2[1], pos2_t2[0]) + 2*np.pi])
    print('angle ships',beta)
    t = np.array([t1, t1, t2, t2])
    t_sat = np.array([ts1_t1, ts2_t1, ts1_t2, ts2_t2])

    distance = c * (t - t_sat) #calculated NOT using relativity
    t_sat_rel = np.sqrt(1-2*mass*G/c**2/planet_radius)/np.sqrt(1-2*mass*G/c**2/radius)/np.sqrt(1-vel**2/c**2)*t_sat
    print('KONST WELL', np.sqrt(1-2*mass*G/c**2/planet_radius)/np.sqrt(1-2*mass*G/c**2/radius))
    print('KONST SPEC', 1/np.sqrt(1-vel**2/c**2))

    distance = c * (t - t_sat_rel) #calculated using relativity

    #print(distance)
    alpha = np.arccos((radius**2 + planet_radius**2 - distance**2)/(2*radius*planet_radius))
    print('angle  diff',alpha)
    angles = alpha*np.array([-1,-1,1,1]) - beta*np.array([-1,-1,-1,-1])
    print('ANGLES phone', angles)
    angle1 = angles[0]/2 + angles[1]/2
    angle2 = angles[2]/2 + angles[3]/2
    print(angle1, 'or', angle1*180/np.pi - 180)
    print(angle2, 'or', angle2*180/np.pi - 180)
    position1 = np.array([planet_radius * np.cos(angle1), planet_radius * np.sin(angle1)])
    print(position1)
    position2 = np.array([planet_radius * np.cos(angle2), planet_radius * np.sin(angle2)])
    print(position2)
    print(nt.norm(position2 - position1), 'meters difference')


if __name__ == '__main__':
    #Part9_Ex3_P1()
    #Part9_Ex3_P2()
    Part9_Ex5()

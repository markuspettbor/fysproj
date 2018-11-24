import variables as vars
import numpy as np
import matplotlib.pyplot as plt

solar_system = vars.solar_system
chosen_planet = 1

#solar_system.part2A_4(chosen_planet, friend_seed=None, increase_height=False, filename1="part2A4_frame1.xml", filename2="part2A4_frame2.xml")
#solar_system.part2A_5(chosen_planet,friend_seed=None, filename1="part2A5_frame1.xml", filename2="part2A5_frame2.xml")

def twin_code():
    a = -0.1/vars.c*vars.year #per year
    L0 = 200 #c
    v0 = 0.99 #c
    tB = L0/v0
    tBm = L0/v0 - L0*v0
    t_turn = tB - v0/a
    print(t_turn)
    tY = np.linspace(tB, t_turn, 6001)
    v = v0 + a*(tY-tB)
    xY = L0 + v0*(tY-tB) + 1/2*a*(tY-tB)**2
    tYm = xY*v +tY
    print(v)
    plt.plot(tY, tYm)
    #plt.axis('equal')
    plt.show()
    print('Stopped at t =', t_turn)

if __name__ == '__main__':
    twin_code()

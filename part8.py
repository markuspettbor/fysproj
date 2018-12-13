import variables as vars
import numpy as np
import matplotlib.pyplot as plt

solar_system = vars.solar_system
chosen_planet = 1

#solar_system.part2A_4(chosen_planet, friend_seed=None, increase_height=False, filename1="part2A4_frame1.xml", filename2="part2A4_frame2.xml")
#solar_system.part2A_5(chosen_planet,friend_seed=None, filename1="part2A5_frame1.xml", filename2="part2A5_frame2.xml")
#solar_system.part2B_5(chosen_planet,friend_seed=None, filename1="part2B5_frame1.xml", filename2="part2B5_frame2.xml")

#solar_system.part2B_1(chosen_planet,friend_seed=None, filename1="part2B1_frame1.xml", filename2="part2B1_frame2.xml", filename3="part2B1_frame3.xml")
#solar_system.part2B_4(chosen_planet,friend_seed=None, increase_height=False, filename1="part2B_4_frame1.xml", filename2="part2B_4_frame2.xml")
#solar_system.part2C_5(number_of_light_signals=30, friend_seed=None, consider_light_travel=False, write_text=True, filename1="part2C_5_frame1.xml", filename2="part2C_5_frame2.xml")
#solar_system.part2C_8(chosen_planet, theta=None, friend_seed=None, increase_height=False, filename="part2C_8.xml")

def twin_code():
    a = -0.1/vars.c*vars.year #per year
    L0 = 200 #c
    v0 = 0.99 #c
    tB = L0/v0
    tBm = L0/v0 - L0*v0
    t_turn = tB - v0/a
    print(t_turn, tB)
    tY = np.linspace(tB, t_turn, 6001)
    v = v0 + a*(tY-tB)
    xY = L0 + v0*(tY-tB) + 1/2*a*(tY-tB)**2
    tYm = tY - xY*v
    print(v)
    plt.plot(tY, tYm, '-k', linewidth  = 0.8)
    plt.xlabel('$T_Y$ [years]', size = 12)
    plt.ylabel('$T_{Y\'}$ [years]', size  = 12)
    #plt.axis('equal')
    plt.show()
    print('Stopped at t =', t_turn)

if __name__ == '__main__':
    twin_code()

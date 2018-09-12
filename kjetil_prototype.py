#compare two areas one where close to aphelion, one where close to perihelion
#INPUTS: positions_vect, velocities_vect, time_vect..
import numpy as np
import variables as var
import part2 as p2

time_vect = p2.all_time[0]
position_vect_x = p2.orbital_x[0]
position_vect_y = p2.orbital_y[0]
velocity_vect_x = p2.velocity_vx[0]
velocity_vect_y = p2.velocity_vy[0]
#print(position_vect)
a = var.a[0]
m_star = var.m_star
m = var.m[0]
orbital_period = p2.kepler3(m_star, m, a)
for n in range(len(time_vect)):
    if time_vect[n] > orbital_period:
        i = n
        break
print(i)
position_vect_x = position_vect_x[0:i] #eller motsatt
position_vect_y = position_vect_y[0:i] #eller motsattnoice = np.random.normal(0, max(vr)/5, len(vr)) #mu sigma length
position = np.array([position_vect_x, position_vect_y])
position = position[:,0:10]
print(position)
#print(position_vect_x)
velocity_vect_x = velocity_vect_x[0:i] #eller motsatt
velocity_vect_y = velocity_vect_y[0:i] #eller motsatt
time_vect = time_vect[0:i]

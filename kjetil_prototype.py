#compare two areas one where close to aphelion, one where close to perihelion
#INPUTS: positions_vect, velocities_vect, time_vect..
import numpy as np
import variables as var
import part2 as p2
time_vect = np.linspace(0,10,1000)
positions_vect = orbitals_x[0]# [0] is starting planet
a = var.a[0]
m_star = var.m_star
m = var.m[0]
#orbital_period = kepler3(m_star, m, a)
i = next(x for x in time_vect if x >= orbital_period)
print(i)
positions_vect = positions_vect[:, 0:i] #eller motsatt
print(positions_vect)
velocity_vect = velocity_vect[:, 0:i] #eller motsatt
time_vect = time_vect[:, 0:i]

#compare two areas one where close to aphelion, one where close to perihelion
#INPUTS: positions_vect, velocities_vect, time_vect..
import numpy as np
import matplotlib.pyplot as plt
import variables as var
import numtools as nt

def kepler3(m1, m2, a):
    return np.sqrt(4*np.pi**2/(var.G*(m1+m2))*a**3)


a = var.a
m_star = var.m_star
m = var.m

orbital_period = kepler3(m_star, m, a)
orbital_period = np.append([7], orbital_period)
print(orbital_period)
t = np.load('saved/saved_orbits/launch_resolution/time_onlysun.npy')
x = np.load('saved/saved_orbits/launch_resolution/pos_onlysun.npy')
v = np.load('saved/saved_orbits/launch_resolution/vel_onlysun.npy')
print(x.shape)
dA = np.zeros(len(x[0]))#, len(x[0][0])])
for n in range(int(len(x[0]))):
    start = 10
    if n == 4:
        start = 2000
    vel = v[:,n,start:int(orbital_period[n]/(t[1]-t[0])*1.1)+start]
    pos = x[:,n,start:int(orbital_period[n]/(t[1]-t[0])*1.1)+start]
    if n > 0:
        print('Planet %i:' %(n))
    else:
        print('Sun:')
    dist = nt.norm(pos)      #xfin distance from cm
    apoapsis = np.argmax(dist)   #xfin max distance [index] from cm
    periapsis = np.argmin(dist)   #xfin min distance [index] from cm
    buel_api = nt.norm(pos[:,apoapsis+1] - pos[:,apoapsis-1])
    area_api = dist[apoapsis]*buel_api/2
    #print('area api  =', area_api)
    #print('dist_api  =', buel_api)
    #print('velo_api  =', nt.norm(vel[:,apoapsis]))
    buel_peri = nt.norm(pos[:,periapsis+1] - pos[:,periapsis-1])
    area_peri = dist[periapsis]*buel_peri/2
    #print('\narea peri =', area_peri)
    #print('dist_peri =', buel_peri)
    #print('velo_peri =', nt.norm(vel[:,periapsis]))
    print('Apoapsis / Periapsis =', area_api/area_peri)
    plt.plot(dist)
    print('')

plt.show()

'''
print(apoapsis[0])
#measure area perihelion
for i in num:
    rad = dist[i][apoapsis[i]]
    dx = nt.norm(x[i][apoapsis[i]+1]- x[i][apoapsis[i]])
    dA[i] = rad*dx/2
print(dA[0])
for i in num:
    rad = dist[i][periapsis[i]]
    dx = abs(dist[i][periapsis[i]+1]- dist[i][periapsis[i]])
    dA[i] = rad*dx/2
print(dA[0])#measure area, dA =
'''

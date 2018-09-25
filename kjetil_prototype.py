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
t = np.load('saved_orbits/10k_20o/time.npy')
x = np.load('saved_orbits/10k_20o/pos.npy')
v = np.load('saved_orbits/10k_20o/vel.npy')
print(x.shape)
dA = np.zeros(len(x[0]))#, len(x[0][0])])
for n in range(int(len(x[0]))):
    start = 10
    if n == 4:
        start = 2000
    vel = v[:,n,start:int(orbital_period[n]/(t[1]-t[0])*1.1)+start]
    print(vel.shape)
    pos = x[:,n,start:int(orbital_period[n]/(t[1]-t[0])*1.1)+start]
    print(pos.shape)
    if n > 0:
        print('Planet %i:' %(n))
    else:
        print('Sun:')
    dist = nt.norm(pos)      #xfin distance from cm
    maks = np.argmax(dist)   #xfin max distance [index] from cm
    mini = np.argmin(dist)   #xfin min distance [index] from cm
    buel_api = nt.norm(pos[:,maks+1] - pos[:,maks-1])
    area_api = dist[maks]*buel_api/2
    print('area api  =', area_api)
    print('dist_api  =', buel_api)
    print('velo_api  =', nt.norm(vel[:,maks]))
    buel_peri = nt.norm(pos[:,mini+1] - pos[:,mini-1])
    area_peri = dist[mini]*buel_peri/2
    print('\narea peri =', area_peri)
    print('dist_peri =', buel_peri)
    print('velo_peri =', nt.norm(vel[:,mini]))
    print('\napi/peri =', area_api/area_peri)
    plt.plot(dist)
    print('')

plt.show()

'''
print(maks[0])
#measure area perihelion
for i in num:
    rad = dist[i][maks[i]]
    dx = nt.norm(x[i][maks[i]+1]- x[i][maks[i]])
    dA[i] = rad*dx/2
print(dA[0])
for i in num:
    rad = dist[i][mini[i]]
    dx = abs(dist[i][mini[i]+1]- dist[i][mini[i]])
    dA[i] = rad*dx/2
print(dA[0])#measure area, dA =
'''

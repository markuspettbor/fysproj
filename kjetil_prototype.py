#compare two areas one where close to aphelion, one where close to perihelion
#INPUTS: positions_vect, velocities_vect, time_vect..
import numpy as np
import matplotlib.pyplot as plt
import variables as var
import part2 as p2
import numtools as nt

t = p2.time2[0]
a = var.a[0]
m_star = var.m_star
m = var.m[0]
orbital_period = p2.kepler3(m_star, m, a)
x = p2.xfin[:,1:]
v = p2.vfin[:,1:]
x = x[:,:,100:-100]
v = v[:,:,100:-100]
#x = np.reshape(x, (len(x[0]),len(x),len(x[0][0])))
#v = np.reshape(v, (len(v[0]),len(v),len(v[0][0])))
print(x.shape)
#find index of xfin planet max
#dist = np.zeros([len(x), len(x[0][0])])
#maks = np.zeros(len(x), dtype = 'int')#, len(x[0][0])])
#mini = np.zeros(len(x), dtype = 'int')#, len(x[0][0])])
dA = np.zeros(len(x))#, len(x[0][0])])
num = np.linspace(0,len(x)-1, len(x), dtype = 'int')
for n in range(len(x[0])):
    vel = v[:,n,:]
    pos = x[:,n,:]
    print(pos.shape)
    print('Planet %i:' %(n+1))
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

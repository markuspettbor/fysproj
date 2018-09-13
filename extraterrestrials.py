import variables as var
import numpy as np
import numtools as nt
import matplotlib.pyplot as plt
import part2 as p2

def radial_velocity_curve():
    vel_sun = np.array([vv.transpose()[0,0], vv.transpose()[1,0]])
    inc = 3/7*np.pi
    vel_pec = 9e3
    #print('shape vel =', vel_sun.shape )
    vel_radial_data = nt.create_radial_velocity(vel_sun, vel_pec, inc)
    plt.plot(vel_radial_data)
    plt.show()


def create_light_curve(pos, rp, rs): #position, time, radius planet, radius sun
    #print('shape pos =', pos.shape )
    x = pos[0] #hope this is right
    y = pos[1]
    #plt.plot(x/max(x), 'm')
    #plt.plot(y/max(y), 'g')
    #print(max(x) - rs)
    #plt.plot(y)
    #plt.show()
    area_planet = np.pi*rp**2
    area = np.zeros(len(x))
    n = np.linspace(0,len(x)-1, len(x), dtype=int)
    for x, y, n in zip(x, y, n):
        if y < 0:
            area[n] = 0
        elif abs(x) >= rs+rp:
            area[n] = 0
        elif abs(x) <= rs-rp:
            area[n] = area_planet
        elif nt.norm([x,rp]) > rs:
            if x < 0: #innenfor venstre
                #print('ONE')
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                + c*np.sqrt(rs**2-c**2) + rs**2*np.arcsin(c/rs) + 2*c*x
                area[n] = areal
            else:   #innenfor høyre
                #print('TWO')
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                + c*np.sqrt(rs**2-c**2) + rs**2*np.arcsin(c/rs) - 2*c*x
                area[n] = areal
        else:   #utenfor venstre
            if x < 0:
                #print('THREE')
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal_inv = +c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                - c*np.sqrt(rs**2-c**2) - rs**2*np.arcsin(c/rs) - 2*c*x
                area[n] = area_planet - areal_inv
            else: #utenfor høyre
                #print('FOUR')
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal_inv = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                - c*np.sqrt(rs**2-c**2) - rs**2*np.arcsin(c/rs) + 2*c*x
                area[n] = area_planet - areal_inv
    #plt.plot(area)
    #plt.show()
    return area

m1 = var.m
m2 = var.m_star
orbit = np.sqrt(var.x0**2 + var.y0**2)
radius_planets = var.radius*1000
radius_sun = var.radius_star*1000
index = np.argmax(m1*m2/(orbit**2))
#print(index)
#rp = np.array([var.x0[index],var.y0[index]])
#print(rp)
radius_heavy = radius_planets[index]

mask = p2.mask

xx = p2.xx * var.AU_tall
vv = p2.vv * var.AU_tall/(365.26*24*60*60)
#bodies = len(p2.mass)
# rr = xx[:,1:]
# print('RR', rr.shape)
# plt.plot(rr[0,0], 'k')
# plt.plot(rr[1,0], 'y')
# plt.show()
# print('shape xx =', xx.shape)
# print('max [0][0]', max(xx[0,0]))
# print('max [0][1]', max(xx[0,1]))
# print('max [1][0]', max(xx[1,0]))
# print('max [1][1]', max(xx[1,1]))
pos_planets = np.copy(xx[:,1:])


#print('shape pos_planets =', pos_planets.shape)

area_sun = np.pi*radius_sun**2
area_covered = np.zeros(len(xx[0,0]))
#flux_relative = np.zeros([1, len(xx[0,0])])
#flux_relative_data = 0*(len(xx[0,0]))
#print(area_covered.shape, 'SHAPE')
n = 0
for mask_index in mask:
    # plt.plot(pos_planets[0,n], 'm')
    # plt.plot(pos_planets[1,n], 'g')
    # plt.show()
    # print('max [0]',max(pos_planets[0,n]))
    # print('max [1]',max(pos_planets[1,n]))
    # print('n=', n)
    # print('shape pos_planet =', pos_planets[:,n].shape)
    area_covered += create_light_curve(pos_planets[:,n], radius_planets[mask_index], radius_sun)
    #plt.plot((area_sun - area_covered)/area_sun)
    n+=1
flux_relative = (area_sun-area_covered)/area_sun
flux_relative_data = nt.noiceify(flux_relative, 0, 0.2)
plt.plot(flux_relative_data, 'r')
plt.plot(flux_relative, 'b')
plt.show()


radial_velocity_curve()

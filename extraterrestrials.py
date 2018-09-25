import variables as var
import numpy as np
import numtools as nt
import matplotlib.pyplot as plt
import part2 as p2

def radial_velocity_curve():
    vel_sun = np.copy(vv[:,0])
    inc = 3/7*np.pi
    vel_pec = 0.420
    vel_radial_data = nt.create_radial_velocity(vel_sun, vel_pec, inc)
    return vel_radial_data

def create_light_curve(pos, rp, rs): #position, time, radius planet, radius sun
    x = pos[0]
    y = pos[1]
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
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                + c*np.sqrt(rs**2-c**2) + rs**2*np.arcsin(c/rs) + 2*c*x
                area[n] = areal
            else:   #innenfor høyre
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                + c*np.sqrt(rs**2-c**2) + rs**2*np.arcsin(c/rs) - 2*c*x
                area[n] = areal
        else:   #utenfor venstre
            if x < 0:
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal_inv = +c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                - c*np.sqrt(rs**2-c**2) - rs**2*np.arcsin(c/rs) - 2*c*x
                area[n] = area_planet - areal_inv
            else: #utenfor høyre
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal_inv = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                - c*np.sqrt(rs**2-c**2) - rs**2*np.arcsin(c/rs) + 2*c*x
                area[n] = area_planet - areal_inv
    return area
m1 = var.m
m2 = var.m_star
orbit = np.sqrt(var.x0**2 + var.y0**2)
radius_planets = var.radius*1000/var.AU_tall
radius_sun = var.radius_star*1000/var.AU_tall
index = np.argmax(m1*m2/(orbit**2))
radius_heavy = radius_planets[index]
mask = p2.mask
xx = p2.xfin# * var.AU_tall
vv = p2.vfin# * var.AU_tall/(365.26*24*60*60)
time = p2.time2
pos_planets = np.copy(xx[:,1:])
area_sun = np.pi*radius_sun**2
area_covered = np.zeros(len(xx[0,0]))
area_covered_saved = []
n = 0
for mask_index in mask:
    area = create_light_curve(pos_planets[:,n], radius_planets[mask_index], radius_sun)
    area_covered += area
    area_covered_saved.append(area)
    n+=1
flux_relative = (area_sun-area_covered)/area_sun
mu = 0
sig = max(area_covered/area_sun*0.2)
flux_relative_data = nt.noiceify(flux_relative, mu, sig)
#start = int(len(flux_relative_data)*(1.681/2))
#stop = int(len(flux_relative_data)*(1.686/2))
#plt.plot(time[start:stop], flux_relative_data[start:stop], 'r')
#plt.plot(time, flux_relative, 'b')
#plt.show()

vel_radial_data = radial_velocity_curve()
plt.plot(time, vel_radial_data)
plt.show()

#def save_data():
    #velocity_radial_data_sending = np.array([time, vel_radial_data]).transpose()
    #np.savetxt('velocity_radial_one_planet_v3.txt', velocity_radial_data_sending)

    #flux_relative_data_sending = np.array([time[start:stop], flux_relative_data[start:stop]]).transpose()
    #np.savetxt('flux_relative_V2.txt', flux_relative_data_sending)
#save_data()

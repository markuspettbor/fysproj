import variables as var
import numpy as np
import numtools as nt
import matplotlib.pyplot as plt
import part2 as p2

m1 = var.m
m2 = var.m_star
orbit = np.sqrt(var.x0**2 + var.y0**2)
radius_temp = var.radius*1000
radius_sun = var.radius_star*1000
index = np.argmax(m1*m2/(orbit**2))
#print(index)
#rp = np.array([var.x0[index],var.y0[index]])
#print(rp)
radius = radius_temp[index]

xx = p2.xx * var.AU_tall #IKKE 110 PROSENT SIKKER
vv = p2.vv * var.AU_tall/(365.26*24*60*60)
pos_plan = (xx[:,1,:][:,0],xx[:,1,:][:,1]) #(xx[:,i,:][:,0],xx[:,i,:][:,1])
vel_sun = (vv[:,0,:][:,0],vv[:,0,:][:,1])

'''
vr1 = np.sin(np.linspace(0, 2*np.pi, 500000))*radius_sun*2
vr2 = np.cos(np.linspace(0, 2*np.pi, 500000))*radius_sun*2
vr  = np.array([vr1,vr2])
tid = np.linspace(0,1,500)
'''
tid = 0 #TEMPORARY
i = 3/7*np.pi
vel_pec = 10e3
vel_radial_data = nt.create_radial_velocity(vel_sun, vel_pec, i)


def create_light_curve(pos, t, rp, rs): #position, time, radius planet, radius sun
    x = pos[1] #hope this is right
    y = pos[0]
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
            if x < 0:
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                + c*np.sqrt(rs**2-c**2) + rs**2*np.arcsin(c/rs) + 2*c*x
                area[n] = areal
            else:
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                + c*np.sqrt(rs**2-c**2) + rs**2*np.arcsin(c/rs) - 2*c*x
                area[n] = areal
        else:
            if x < 0:
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal_inv = +c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                - c*np.sqrt(rs**2-c**2) - rs**2*np.arcsin(c/rs) - 2*c*x
                area[n] = area_planet - areal_inv
            else:
                c = np.sqrt(rs**2-((rs**2-rp**2+x**2)/(2*x))**2) #c = y i kryssningspunkt mellom sirklene
                areal_inv = c*np.sqrt(rp**2-c**2) + rp**2*np.arcsin(c/rp) \
                - c*np.sqrt(rs**2-c**2) - rs**2*np.arcsin(c/rs) + 2*c*x
                area[n] = area_planet - areal_inv
    #plt.plot(area)
    #plt.show()
    return area
area_covered = create_light_curve(pos_plan, tid, radius, radius_sun)
area_sun = np.pi*radius_sun**2
flux_relative = (area_sun-area_covered)/area_sun
flux_relative_data = nt.noiceify(flux_relative, 0, 0.2)
plt.plot(vel_radial_data)
plt.show()
plt.plot(flux_relative_data, 'r')
plt.plot(flux_relative, 'b')
plt.show()

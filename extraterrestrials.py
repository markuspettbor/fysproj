import variables as var
import numpy as np
import numtools as nt
import matplotlib.pyplot as plt
import part2 as p2

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
bodies = len(p2.mass)
pos_plan = np.zeros(bodies)
for i in range(bodies):
    pos_plan[i] = (xx[0, i+1], xx[1,i+1]) #excludes sun with + 1

area_sun = np.pi*radius_sun**2
flux_relative = np.zeros(bodies-1)
flux_relative_data = np.zeros(bodies-1)

for mask_index, pos_plan in zip(mask, pos_plan)
    area_covered = create_light_curve(pos_plan, radius[mask_index], radius_sun)
    flux_relative[i] = (area_sun-area_covered)/area_sun
    flux_relative_data[i] = nt.noiceify(flux_relative, 0, 0.2)
plt.plot(flux_relative_data[0], 'r')
plt.plot(flux_relative[0], 'b')
plt.show()


radial_velocity_curve():

def radial_velocity_curve():
    vel_sun = (vv[0,0], vv[1,0])
    inc = 3/7*np.pi
    vel_pec = 9e3
    vel_radial_data = nt.create_radial_velocity(vel_sun, vel_pec, inc)
    plt.plot(vel_radial_data)
    plt.show()


def create_light_curve(pos, rp, rs): #position, time, radius planet, radius sun
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

import numpy as np
import orbit_tools as ot
import numtools as nt
import variables as vars
import launch
import matplotlib.pyplot as plt

def solar_panel_area(watt, effic, F):
    return watt/(F*effic)*np.pi/2
radius_star = vars.radius_star*1000
dist_AU = np.sqrt(vars.x0**2 + vars.y0**2)
sigma = vars.sbc
T = vars.temp
radius = vars.radius_normal_unit #radius of planet
dist = dist_AU*vars.AU_tall
F = sigma*T**4 #flux per square meter of sun
surface_sun = 4*np.pi*(radius_star)**2
surface_shells = 4*np.pi*dist**2
L_tot = F*surface_sun #total W from the sun

F_rec = L_tot/surface_shells #W/m²
F_rec = F*surface_sun/surface_shells
solar_panel_A = solar_panel_area(40, 0.12, F_rec)
planet_temperature = T*(radius_star**2/dist**2/4)**(1/4)
planet_temperature = (F_rec/4/sigma)**(1/4)
#planet_temperature = (1/4*T**4*surface_sun/surface_shells)**(1/4)
Temps = np.array([390, 260])
distanser = np.sqrt(T**4/Temps**4*radius_star**2/4)
print(distanser/vars.AU_tall)
print(dist)


if __name__ == '__main__':
    #print('At a distance %.6f AU from the sun, the flux equals %.3f W/m², and the lander needs %.3f m² of solar panels to function' %(dist_AU[0], F_rec, solar_panel_A))
    for i in [4,0,1,6,5,3,2]:
        print('\nPlanet number %i' %i)
        print('Solar panel area  ', solar_panel_A[i])
        print('Planet temperature', planet_temperature[i]-273.15)

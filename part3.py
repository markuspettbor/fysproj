import numpy as np

import variables as vars

def solar_panel_area(watt, effic):
    return watt/(F_rec*effic)

def temp_calc(i, L):
    r = np.sqrt(vars.x0[i]**2 + vars.y0[i]**2)*vars.AU_tall
    radius = vars.radius_normal_unit[i] #radius of planet
    F_rec = L/(4*np.pi*(r)**2) #W/m²
    energy_recieved = F_rec*np.pi*radius**2
    area = 4*np.pi*radius**2
    return (energy_recieved/(sigma*area))**(1/4)

dist_AU = np.sqrt(vars.x0[0]**2 + vars.y0[0]**2)
T = vars.temp
sigma = vars.sbc
radius = vars.radius_normal_unit #radius of planet
r = dist_AU*vars.AU_tall
F = sigma*T**4 #flux per square meter of sun
surface_sun = 4*np.pi*(vars.radius_star*1000)**2
L_tot = F*surface_sun #total W from the sun
F_rec = L_tot/(4*np.pi*(r)**2) #W/m²
solar_panel_A = solar_panel_area(40, 0.12)
print('At a distance %.6f AU from the sun, the flux equals %.3f W/m², and the lander needs %.3f m² of solar panels to function' %(dist_AU, F_rec, solar_panel_A))

planet_temperature = np.zeros(vars.n)
for i in [4,0,1,6,5,3,2]:
    print('\nPlanet number %i' %i)
    solar_panel_A = solar_panel_area(40, 0.12)
    planet_temperature[i] = temp_calc(i, L_tot)
    print(solar_panel_A)
    print(planet_temperature[i]-273.15)

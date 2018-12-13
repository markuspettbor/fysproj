from ast2000solarsystem import AST2000SolarSystem
import sys, os
import numpy as np

user = 'markusbpkjetilmg'
mmH2 = 2.016 #g/mol
mol = 6.022140857e23 #1/mol
'''https://physics.nist.gov/cgi-bin/cuu/Value?na'''
k = 1.38064852e-23 # Boltzmann constant
molecule_mass = mmH2/mol/1000 #mass of single molecule
gravitational_constant = 6.67408e-11

sys.stdout = open(os.devnull, 'w')
seed = AST2000SolarSystem.get_seed(user)
solar_system = AST2000SolarSystem(seed)
ref_stars = solar_system.get_ref_stars()
sys.stdout = sys.__stdout__

solar_system = AST2000SolarSystem(seed)
n = solar_system.number_of_planets
x0 = solar_system.x0
y0 = solar_system.y0
vx0 = solar_system.vx0
vy0 = solar_system.vy0
a = solar_system.a
e = solar_system.e
theta0 = solar_system.omega
psi0 = solar_system.psi
rho0 = solar_system.rho0
radius = solar_system.radius #radius of planets as list/array
m_star = solar_system.star_mass
radius_star = solar_system.star_radius
m = solar_system.mass #mass of planets as list/array
satellite = solar_system.mass_sat
period = solar_system.period
temp = solar_system.temperature
mass_lander = solar_system.mass_lander
#area_lander = solar_system.area_lander #LIE
area_lander = 0.3
area_sat = solar_system.area_sat
G = 4*np.pi**2
G_SI = 6.67408*1e-11
c = 299792458 #speed of light
'''https://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=speed+of+light'''
solmasse = 1.989e30
AU_tall = 149597871000 #m/AU
m_normal_unit = m*solmasse
radius_normal_unit = radius*1000
radius_AU = radius_normal_unit / AU_tall
year = 365.26*24*60*60
sbc = 5.670367e-8 #Stefan-Boltzmann constant
solar_mass = 1.989e30
'''https://physics.nist.gov/cgi-bin/cuu/Value?sigma'''


if __name__ =='__main__':
    def kepler3(m1, m2, a):
        return np.sqrt(4*np.pi**2/(G*(m1+m2))*a**3)
    for i in [4,0,1,6,5,3,2]:
        print('Planet number %i' %i)
        print('Radius %f' %(radius_normal_unit[i]))
        print('Planet distance %f' %(np.sqrt(x0[i]**2 + y0[i]**2)*AU_tall))
        print('Planet mass %e' %(m[i]*solmasse))
        print('Atmosphere density %f' %(rho0[i]))
        print('Days %f' %(period[i]))
        print('Year %f' %(kepler3(m_star, m[i], a[i])))
        print('Gravity %f' %(G*m_star*m[i]/radius_AU[i]**2/AU_tall*year))
        print('Density %f' %((m_normal_unit[i]/(4/3*np.pi*radius_normal_unit[i]**3))))
        #print('exitrentricity %f' %(e[i]))
        print('')
    print('\n')

    for i in [4,0,1,6,5,3,2]:
        print('Planet number %i' %i)
        print('Radius (AU)%f' %(radius_normal_unit[i]/AU_tall))
        print('Planet distance (AU)%f' %(np.sqrt(x0[i]**2 + y0[i]**2)))
        print('Planet mass (SM)%e' %(m[i]))
        print('Atmosphere density %f' %(rho0[i]))
        print('Days %f' %(period[i]))
        print('Year %f' %(kepler3(m_star, m[i], a[i])))
        print('Gravity %f' %(G*m_star*m[i]/radius_AU[i]**2/AU_tall*year))
        print('Density %f' %((m_normal_unit[i]/(4/3*np.pi*radius_normal_unit[i]**3))))
        #print('exitrentricity %f' %(e[i]))
        print('')
    #print('SEED', seed)

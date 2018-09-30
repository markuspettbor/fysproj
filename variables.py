from AST2000SolarSystem import AST2000SolarSystem
import sys, os
import numpy as np

user = 'markusbpkjetilmg'
mmH2 = 2.016 #g/mol
mol = 6.022140857e23 #1/mol
k = 1.38064852e-23 # Boltzmann constant
molecule_mass = mmH2/mol/1000 #mass of single molecule
gravitational_constant = 6.67408e-11

sys.stdout = open(os.devnull, 'w')
seed = AST2000SolarSystem.get_seed(user)
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
G = 4*np.pi**2
AU_tall = 149597871000 #m/AU
m_normal_unit = m*1.989e30
radius_normal_unit = radius*1000
radius_AU = radius_normal_unit / AU_tall
year = 365*24*60*60
sbc = 5.670367e-8 #Stefan-Boltzmann constant
solar_mass = 1.989e30 
'''https://physics.nist.gov/cgi-bin/cuu/Value?sigma'''

if __name__ =='__main__':
    def kepler3(m1, m2, a):
        return np.sqrt(4*np.pi**2/(G*(m1+m2))*a**3)

    for i in [4,0,1,6,5,3,2]:
        print('Planet number %i' %i)
        print('Planet distance/dumstance %f' %(np.sqrt(x0[i]**2 + y0[i]**2)/np.sqrt(x0[0]**2 + y0[0]**2)))
        print('Planet mass/dummass %f' %(m[i]/m[0]))
        print('Atmosphere density/atmosdumsity %f' %(rho0[i]/rho0[0]))
        print('Days/dumys %f' %(period[i]/period[0]))
        print('Year/dumears %f' %(kepler3(m_star, m[i], a[i])/kepler3(m_star, m[0], a[0])))
        print('Gravity/dumvity %f' %(G*m_star*m[i]/radius_AU[i]**2/AU_tall*year/(G*m_star*m[0]/radius_AU[0]**2/AU_tall*year)))
        print('density/dumsity %f' %((m_normal_unit[i]/(4/3*np.pi*radius_normal_unit[i]**3))/(m_normal_unit[0]/(4/3*np.pi*radius_normal_unit[0]**3))))
        #print('exitrentricity %f' %(e[i]))
        print('')

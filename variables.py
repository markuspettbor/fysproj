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
radius = solar_system.radius #radius of planets as list/array
m_star = solar_system.star_mass
m = solar_system.mass #mass of planets as list/array
G = 4*np.pi**2
m_normal_unit = m*1.989e30
radius_normal_unit = radius*1000

satelite = solar_system.mass_sat

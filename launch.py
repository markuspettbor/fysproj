import variables
import numpy as np

radius = variables.radius_normal_unit[0]
planet_mass = variables.m_normal_unit[0]
print(planet_mass)
print(radius)
grav_const = variables.gravitational_constant
position = radius
force = 813e3
fuel_consumption = force * 1e-4
initial_mass = 200e3
mass = initial_mass
dt = 0.001
t = 0
escape_velocity = (2*grav_const*planet_mass/position)**0.5
print(escape_velocity)
velocity = 0
while velocity < escape_velocity: #1Dimentional
    acceleration = force/mass# - grav_const*(planet_mass*mass)/position**2
    velocity = velocity + acceleration*dt
    position = position + velocity*dt
    mass = mass - fuel_consumption*dt
    print(position)
    escape_velocity = np.sqrt(2*grav_const*planet_mass/position)
    t = t + dt
print(t/60)

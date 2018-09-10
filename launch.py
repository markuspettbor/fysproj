import matplotlib.pyplot as plt
import variables
def launch(force_box, fuel_box):
    radius = variables.radius_normal_unit[0]
    planet_mass = variables.m_normal_unit[0]
    grav_const = variables.gravitational_constant
    satellite_mass = variables.satellite
    position = radius
    rocket_mass = 0 #input?
    force = 1680e3 #input? regne ut?
    fuel_mass = 100e3 #input? regne ut?
    dt = 0.01 #input?
    t0 = 0 #input?
    boxes = force/force_box
    fuel_consumption = boxes * fuel_box
    mass = satellite_mass + fuel_mass + rocket_mass
    initial_mass = mass
    initial_fuel_mass = fuel_mass
    t = t0
    escape_velocity = (2*grav_const*planet_mass/position)**0.5
    velocity = 0; count = 0; has_fuel = 1
    print('Escape Velocity at Surface = %.3e' % escape_velocity)
    print('Planet Mass = %.3e' % planet_mass)
    print('Plannet Radius = %.3e' % radius)
    print('g =',grav_const*planet_mass/(position**2))
    print('Initial Mass of Rocket = %.3e' % initial_mass)
    while velocity < escape_velocity: #1Dimentional
        acceleration = force/mass*has_fuel - grav_const*planet_mass/(position**2)
        #print(acceleration)
        velocity = velocity + acceleration*dt
        position = position + velocity*dt
        mass = mass - fuel_consumption*dt
        fuel_mass = fuel_mass - fuel_consumption*dt
        #print(position)
        if count % 100 == 0:
            plt.scatter(t, acceleration)
        escape_velocity = (2*grav_const*planet_mass/position)**0.5
        if fuel_mass < 0:
            has_fuel = 0
        t = t + dt
        count += 1
    print('Final Mass of Rocket = %.3e' % mass)
    print('Fuel Percentage Left = %.3f' % (100*fuel_mass/initial_fuel_mass))
    print('Launch Time = %.3f Minutes' % (t/60))
    print('Final Position = %.3e' % position)
    print('Final Velocity = %.3e' % velocity)
    print('Boxes Used = %.3e' % boxes)
    plt.show()

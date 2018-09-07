import matplotlib.pyplot as plt
import variables
def launch(force_box, fuel_box):
    radius = variables.radius_normal_unit[0]
    planet_mass = variables.m_normal_unit[0]
    print('Planet Mass =', planet_mass)
    print('Plannet Radius =', radius)
    grav_const = variables.gravitational_constant
    position = radius
    force = 1680e3
    boxes = force/force_box
    fuel_consumption = boxes * fuel_box
    initial_mass = 200e3
    mass = initial_mass
    print('g =',grav_const*planet_mass/(position**2))
    dt = 0.01
    t = 0
    escape_velocity = (2*grav_const*planet_mass/position)**0.5
    print('Escape Velocity at Surface =', escape_velocity)
    velocity = 0
    count = 0
    print('Initial Mass of Rocket =' , initial_mass)
    while velocity < escape_velocity: #1Dimentional
        acceleration = force/mass - grav_const*planet_mass/(position**2)
        #print(acceleration)
        velocity = velocity + acceleration*dt
        position = position + velocity*dt
        mass = mass - fuel_consumption*dt
        #print(position)
        if count % 100 == 0:
            plt.scatter(t, acceleration)
        escape_velocity = (2*grav_const*planet_mass/position)**0.5
        t = t + dt
        count = count+1
    print('Final Mass of Rocket =', mass)
    print('Fuel Percentage Left =', mass/initial_mass)
    print(t/60)
    print('Final Position =', position)
    print('Final Velocity =', velocity)
    print('Boxes Used =', boxes)
    plt.show()

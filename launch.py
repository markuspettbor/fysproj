import matplotlib.pyplot as plt
import numtools as nt
import numpy as np
import variables
def launch(force_box, fuel_box, testing = False):
    radius = variables.radius_normal_unit[0]
    planet_mass = variables.m_normal_unit[0]
    grav_const = variables.gravitational_constant
    satellite_mass = variables.satellite
    period = variables.period[0] #unit is 24h
    position = radius
    rocket_mass = 0 #input?
    force = 1000e3 #input? regne ut?
    fuel_mass = 80e3 #input? regne ut?
    dt = 0.01 #input?
    t0 = 0 #input?
    phi = 0*np.pi/(period*24*60*60)
    theta = 0*np.pi
    boxes = force/force_box
    fuel_consumption = boxes * fuel_box
    mass = satellite_mass + fuel_mass + rocket_mass
    initial_mass = mass
    initial_fuel_mass = fuel_mass
    t = t0
    escape_velocity = (2*grav_const*planet_mass/position)**0.5
    rot_velocity = 0 #2*np.pi*radius/(period*24*60*60)
    velocity = 0; count = 0; has_fuel = 1
    part_consumed_box = fuel_box/variables.molecule_mass
    print('Escape Velocity at Surface = %.3e' % escape_velocity)
    print('Planet Mass = %.3e' % planet_mass)
    print('Plannet Radius = %.3e' % radius)
    print('g =',grav_const*planet_mass/(position**2))
    print('Initial Mass of Rocket = %.3e' % initial_mass)
    while (velocity**2 + rot_velocity**2)**0.5 < escape_velocity: #1Dimentional
        acceleration = force/mass*has_fuel - grav_const*planet_mass/(position**2)
        velocity = velocity + acceleration*dt
        position = position + velocity*dt
        mass = mass - fuel_consumption*dt
        fuel_mass = fuel_mass - fuel_consumption*dt
        if count % 100 == 0:
            plt.scatter(t, acceleration)
        escape_velocity = (2*grav_const*planet_mass/position)**0.5
        if fuel_mass < 0:
            has_fuel = 0
            print('NO ORBIT')
            break
        t = t + dt
        count += 1
    velocity = (velocity**2 + rot_velocity**2)**0.5
    print('Final Mass of Rocket = %.3e' % mass)
    print('Fuel Percentage Left = %.3f' % (100*fuel_mass/initial_fuel_mass))
    print('Launch Time = %.3f Minutes' % (t/60))
    print('Final Position = %.3e' % position)
    print('Final Velocity = %.3e' % velocity)
    print('Boxes Used = %.3e' % boxes)
    plt.show()

    #solar system frame of reference
    #rotasjon om planeten med vinkel phi = 2*np.pi/(period*24*60*60)
    #launchsite theta = vinkel til launchsite
    #position_planet, velocity_planet = orbit_calculation(t0, t) #where t0 is start of launch and t is launch duration
        #returnerer ny position til planet om sola, og hastighet om sola
    print([position, 0])
    phi = t*np.pi/(period*24*60*60)
    position_rs = nt.rotate(np.array([position,0]).transpose(), phi + theta) #returns array with [x, y]
    velocity_rs = nt.rotate(np.array([velocity,0]).transpose(), phi + theta)
    print(position_rs)
    print(velocity_rs)
    position_planet, velocity_planet = np.array([variables.x0[0], variables.y0[0]]), np.array([variables.vx0[0], variables.vy0[0]]) #where t0 is start of launch and t is launch duration
    position_final = position_planet*variables.AU_tall + position_rs # position_planet given from orbital calculation
    velocity_final = velocity_planet + velocity_rs # position_planet given from orbital calculation
    #assume all previous is given in SI units
    position_final_AU = position_final/variables.AU_tall #AU
    velocity_final_AU = velocity_final*60*60*24*365/variables.AU_tall #AU/year
    time_AU = t/(60*60*24*365)

    if testing == True:
        variables.solar_system.engine_settings(force_box, boxes, part_consumed_box, initial_fuel_mass, \
        t, nt.rotate(np.array([variables.x0[0] + variables.radius_AU[0], 0]).transpose(), theta), 0)
        variables.solar_system.mass_needed_launch(position_final_AU , test = True)

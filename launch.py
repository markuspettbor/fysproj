import matplotlib.pyplot as plt
import numtools as nt
import numpy as np
import variables as vars
import part4_kjetil as p4k
from scipy.interpolate import interp1d

def launch(time_vector, planet_position, planet_velocity, t0, theta = -1/2*np.pi, testing = False):
    #theta = launchsite on planet, might want to correct for phi rotation before setting theta.
    radius = vars.radius_normal_unit[0]
    planet_mass = vars.m_normal_unit[0]
    grav_const = vars.gravitational_constant
    satellite_mass = vars.satellite
    period = vars.period[0] #unit is 24h
    force_box = np.load('saved/engine/force_box.npy')
    fuel_consumed_box_per_sec = np.load('saved/engine/fuel_consumed_box_per_sec.npy')
    part_consumed_box = fuel_consumed_box_per_sec/vars.molecule_mass
    position = radius
    x_interp = nt.interp_xin(time_vector, planet_position[:,1])
    v_interp = nt.interp_xin(time_vector, planet_velocity[:,1])

    force = 35000e3 #input? kraft for hele raketten
    fuel_mass = 3000e3-satellite_mass #input? total fuelmasse for raketten
    boxes = force/force_box #antall bokser regnes utifra ønsket kraft
    fuel_consumption = boxes * fuel_consumed_box_per_sec #drivstofforbruk regnes ut for samlet motor
    mass = satellite_mass + fuel_mass
    initial_mass = mass
    initial_fuel_mass = fuel_mass
    dt = 0.01
    t = t0
    escape_velocity = (2*grav_const*planet_mass/position)**0.5
    velocity = 0; count = 0; has_fuel = 1

    while velocity < escape_velocity: #1Dimentional
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

    position = position/vars.AU_tall
    velocity = velocity/vars.AU_tall*vars.year
    t_AU = t/vars.year
    phi = t*np.pi*2/(period*24*60*60) #rotasjon til planet om seg selv
    planet_position_t1 = np.array([x_interp[0](t_AU), x_interp[1](t_AU)])
    planet_velocity_t1 = np.array([v_interp[0](t_AU), v_interp[1](t_AU)])
    position_before_rotation = position * nt.unit_vector(planet_position_t1) #radially outwards from sun/cm
    velocity_before_rotation = velocity * nt.unit_vector(planet_position_t1) #radially outwards from sun/cm
    position_after_rotation = nt.rotate(position_before_rotation, phi + theta) #returns array with [x, y]
    velocity_after_rotation = nt.rotate(velocity_before_rotation, phi + theta)
    final_position = planet_position_t1 + position_after_rotation
    final_velocity = planet_velocity_t1 + velocity_after_rotation

    if testing == True:
        vars.solar_system.engine_settings(force_box, boxes, part_consumed_box, initial_fuel_mass, \
        t, nt.rotate(np.array([vars.x0[0] + vars.radius_AU[0], 0]).transpose(), theta), t0/vars.year)

        vars.solar_system.mass_needed_launch(final_position, test = True)

    return t_AU, final_position, final_velocity, mass/vars.solmasse, fuel_mass/vars.solmasse
    #vel and pos relative to planet. Add planets pos and vel after t_AU years to returned values

def test(testing = False):
    t_load = np.load('saved/saved_orbits/launch_resolution/time.npy')
    x_load = np.load('saved/saved_orbits/launch_resolution/pos.npy')
    v_load = np.load('saved/saved_orbits/launch_resolution/vel.npy')
    dt = t_load[1] - t_load[0]

    #print(planet_pos_t0)
    t, pos, vel, mass, fuel_mass = launch(t_load, x_load, v_load, 0, 0, True)
    t1_index = int(t/dt)
    print(t*vars.year)
    print(pos*vars.AU_tall)
    print(vars.radius[0]*1000)
    print(vel*vars.AU_tall/vars.year)
    print(mass*vars.solmasse)
    print(fuel_mass*vars.solmasse)

    if testing == True:
        #index = min(range(len(t_load)), key=lambda i: abs(t_load[i]-t_AU))
        print('INDEX', t1_index)
        x = p4k.position_from_objects(t1_index, vars.solar_system.analyse_distances(), x_load)
        print(x, 'position after launch from part4')

        v = p4k.velocity_from_stars(vars.solar_system.measure_doppler_shifts())
        print(v, 'velocity after launch from part4')
        print(vel, 'velocity from launch.py')




if __name__ == '__main__':
    #force_box = np.load('saved/engine/force_box.npy')
    #fuel_consumed_box_per_sec = np.load('saved/engine/fuel_consumed_box_per_sec.npy')
    #launch(force_box, fuel_consumed_box_per_sec, testing = True)
    test(True)
    pass

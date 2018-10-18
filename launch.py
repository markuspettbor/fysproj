import matplotlib.pyplot as plt
import numtools as nt
import numpy as np
import variables as vars
import part4_kjetil as p4k
from scipy.interpolate import interp1d

def launch(time_vector, planet_position, planet_velocity, t0, theta = -1/2*np.pi, testing = False):
    radius = vars.radius_normal_unit[0]
    planet_mass = vars.m_normal_unit[0]
    grav_const = vars.gravitational_constant
    satellite_mass = vars.satellite
    period = vars.period[0] #unit is 24h
    force_box = np.load('saved/engine/force_box.npy')
    fuel_consumed_box_per_sec = np.load('saved/engine/fuel_consumed_box_per_sec.npy')
    part_consumed_box = fuel_consumed_box_per_sec/vars.molecule_mass
    position = radius
    DDt = time_vector[1]- time_vector[0]
    indd = int(t0/DDt)
    home_planet_position = planet_position[:,1]# - planet_position[:,0]
    home_planet_velocity = planet_velocity[:,1]# - planet_velocity[:,0]

    x_interp = nt.interp_xin(time_vector, home_planet_position)
    v_interp = nt.interp_xin(time_vector, home_planet_velocity)
    planet_position_t0 = np.array([x_interp[0](t0), x_interp[1](t0)])
    planet_velocity_t0 = np.array([v_interp[0](t0), v_interp[1](t0)])

    #print(planet_position_t0)
    #planet_position_t0 = np.array([home_planet_position[0,indd], home_planet_position[1,indd]])
    #planet_velocity_t0 = np.array([home_planet_velocity[0,indd], home_planet_velocity[1,indd]])
    #print(planet_position_t0)

    force = 35000e3 #input? kraft for hele raketten
    fuel_mass = 3000e3-satellite_mass #input? total fuelmasse for raketten
    boxes = force/force_box #antall bokser regnes utifra ønsket kraft
    fuel_consumption = boxes * fuel_consumed_box_per_sec #drivstofforbruk regnes ut for samlet motor
    mass = satellite_mass + fuel_mass
    initial_mass = mass
    initial_fuel_mass = fuel_mass
    dt = 0.01
    t = t0*vars.year
    escape_velocity = (2*grav_const*planet_mass/position)**0.5
    velocity = 0; count = 0; has_fuel = 1
    angular_velocity = 2*np.pi*radius/(period*24*60*60)
    while (velocity**2 + angular_velocity**2)**(1/2) < escape_velocity: #1Dimentional
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
    angular_velocity = angular_velocity/vars.AU_tall*vars.year
    position = position/vars.AU_tall
    velocity = velocity/vars.AU_tall*vars.year
    t_AU = t/vars.year
    planet_position_t1 = np.array([x_interp[0](t_AU), x_interp[1](t_AU)])
    planet_velocity_t1 = np.array([v_interp[0](t_AU), v_interp[1](t_AU)])
    unitvector = nt.unit_vector(planet_position_t1)
    position_before_rotation = position * unitvector + nt.rotate(unitvector, np.pi/2) * (t_AU - t0) * angular_velocity #radiallyish outwards from sun/cm
    velocity_before_rotation = velocity * unitvector + nt.rotate(unitvector, np.pi/2) * angular_velocity
    position_after_rotation = nt.rotate(position_before_rotation, theta) #returns array with [x, y]
    velocity_after_rotation = nt.rotate(velocity_before_rotation, theta)
    final_position = planet_position_t1 + position_after_rotation
    final_velocity = planet_velocity_t1 + velocity_after_rotation
    phi = nt.angle_between(planet_velocity_t1, velocity_after_rotation)
    print('PHI', phi)
    print('THETA', theta)
    if testing == True:
        vars.solar_system.engine_settings(force_box, boxes, part_consumed_box, initial_fuel_mass, \
        t-t0*vars.year, planet_position_t0 + nt.rotate(vars.radius_AU[0]*unitvector, theta), t0)

        vars.solar_system.mass_needed_launch(final_position, test = True)

    return t_AU, final_position, final_velocity, mass/vars.solmasse, fuel_mass/vars.solmasse, phi
    #vel and pos relative to planet. Add planets pos and vel after t_AU years to returned values

def test(testing = False):
    t_load = np.load('saved/saved_orbits/launch_resolution/time_onlysun.npy')
    print(t_load)
    x_load = np.load('saved/saved_orbits/launch_resolution/pos_onlysun.npy')
    v_load = np.load('saved/saved_orbits/launch_resolution/vel_onlysun.npy')
    dt = t_load[1] - t_load[0]
    #print(planet_pos_t0)
    t0 = 0
    t, pos, vel, mass, fuel_mass, phi = launch(t_load, x_load, v_load, t0, -1/2*np.pi, False)
    t, pos, vel, mass, fuel_mass, phi = launch(t_load, x_load, v_load, t0, -1/2*np.pi-phi, True)
    t1_index = int((t)/dt)
    # print(t*vars.year)
    # print(pos*vars.AU_tall)
    # print(vars.radius[0]*1000)
    # print(vel*vars.AU_tall/vars.year)
    # print(mass*vars.solmasse)
    # print(fuel_mass*vars.solmasse)
    if testing == True:
        #index = min(range(len(t_load)), key=lambda i: abs(t_load[i]-t_AU))
        print('INDEX', t1_index)
        x = p4k.position_from_objects(t1_index, vars.solar_system.analyse_distances(), x_load)
        print(x, 'position after launch from part4')
        print(pos, 'position from launch')

        v = p4k.velocity_from_stars(vars.solar_system.measure_doppler_shifts())
        print(v, 'velocity after launch from part4')
        print(vel, 'velocity from launch.py')
        print(nt.norm(x-pos))
        print(nt.norm(v-vel))




if __name__ == '__main__':
    #force_box = np.load('saved/engine/force_box.npy')
    #fuel_consumed_box_per_sec = np.load('saved/engine/fuel_consumed_box_per_sec.npy')
    #launch(force_box, fuel_consumed_box_per_sec, testing = True)
    test(True)
    pass

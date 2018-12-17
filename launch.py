import matplotlib.pyplot as plt
import numtools as nt
import numpy as np
import variables as vars
import part4 as p4
from scipy.interpolate import interp1d
import rocket

def boost(initial_mass, delta_speed_desired, dt = 0.0001):
    satellite_mass = vars.satellite
    thrust = 35000e3 #input? kraft for hele raketten
    fuel_mass = 3000e3-satellite_mass #input? total fuelmasse for raketten
    force_box = np.load('saved/engine/force_box.npy')
    fuel_consumed_box_per_sec = np.load('saved/engine/fuel_consumed_box_per_sec.npy')
    part_consumed_box = fuel_consumed_box_per_sec/vars.molecule_mass
    boxes = thrust/force_box #antall bokser regnes utifra ønsket kraft
    consumption = boxes * fuel_consumed_box_per_sec #drivstofforbruk regnes ut for samlet motor
    mass = initial_mass
    delta_speed = 0
    while delta_speed < delta_speed_desired and mass > satellite_mass:
        delta_speed, mass = nt.euler_fuel_consumption(delta_speed, mass, thrust, consumption, dt)
        #print(delta_speed)
    if mass <= satellite_mass:
        print('out of fuel')
    return mass, delta_speed

def get_engine_settings(t_launch, t_finished):
    '''
    Returns requered launch parameters for engine, except for position and
    launch time (assumed to be known already)
    '''
    satellite_mass = vars.satellite
    force_box = np.load('saved/engine/force_box.npy')
    fuel_consumed_box_per_sec = np.load('saved/engine/fuel_consumed_box_per_sec.npy')
    n_per_box_s = fuel_consumed_box_per_sec/vars.molecule_mass
    force = 35000e3 #input? kraft for hele raketten
    fuel_mass = 3000e3-satellite_mass #input? total fuelmasse for raketten
    boxes = force/force_box #antall bokser regnes utifra ønsket kraft
    fuel_consumption = boxes * fuel_consumed_box_per_sec #drivstofforbruk regnes ut for samlet motor
    mass = satellite_mass + fuel_mass
    initial_mass = mass
    launch_dur = (t_finished - t_launch)*vars.year
    return force_box, boxes, n_per_box_s, fuel_mass, launch_dur,

def launch(time_vector, planet_position, planet_velocity, t0, theta = 1/2*np.pi, testing = False):
    #print('angle = ', theta)
    #print('angle = ', theta * 180/np.pi)
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
    home_planet_position = planet_position[:,1] - planet_position[:,0]
    home_planet_velocity = planet_velocity[:,1] - planet_velocity[:,0]

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
    #print('PARTICLES PER BOX PER SEC', part_consumed_box)
    #print('FUEL CONSUMED PER SEC', fuel_consumption)
    #print('BOXES', boxes)
    dt = 0.01
    t = t0*vars.year
    escape_velocity = (2*grav_const*planet_mass/position)**0.5
    velocity = 0; count = 0; has_fuel = 1
    tangential_velocity = 2*np.pi*radius/(period*24*60*60)
    while (velocity**2 + tangential_velocity**2)**(1/2) < escape_velocity: #1Dimensional
        acceleration = force/mass*has_fuel - grav_const*planet_mass/(position**2)
        velocity = velocity + acceleration*dt
        position = position + velocity*dt
        mass = mass - fuel_consumption*dt
        fuel_mass = fuel_mass - fuel_consumption*dt
        escape_velocity = (2*grav_const*planet_mass/position)**0.5
        if fuel_mass < 0:
            has_fuel = 0
            print('NO ORBIT')
            break
        t = t + dt
        count += 1

    print('\n\n')
    print('Velocity in planet reference frame', velocity)
    print('Duration', t)
    print('Fuel percentage', fuel_mass / initial_fuel_mass*100)
    print('Initial Fuel', initial_fuel_mass)
    print('Fuel', fuel_mass)
    print('height', position-radius)
    print('\n\n')

    tangential_velocity = tangential_velocity/vars.AU_tall*vars.year
    position = position/vars.AU_tall
    velocity = velocity/vars.AU_tall*vars.year
    t_AU = t/vars.year
    planet_position_t1 = np.array([x_interp[0](t_AU), x_interp[1](t_AU)])
    planet_velocity_t1 = np.array([v_interp[0](t_AU), v_interp[1](t_AU)])
    unitvector = nt.unit_vector(planet_position_t1)
    position_before_rotation = position * unitvector + nt.rotate(unitvector, np.pi/2) * (t_AU - t0) * tangential_velocity #radiallyish outwards from sun/cm
    velocity_before_rotation = velocity * unitvector + nt.rotate(unitvector, np.pi/2) * tangential_velocity
    position_after_rotation = nt.rotate(position_before_rotation, theta) #returns array with [x, y]
    velocity_after_rotation = nt.rotate(velocity_before_rotation, theta)
    final_position = planet_position_t1 + position_after_rotation
    final_velocity = planet_velocity_t1 + velocity_after_rotation
    phi = nt.angle_between(planet_velocity_t1, velocity_after_rotation)
    print('PHI', phi)
    print('THETA', theta)
    print('Duration:', t-t0*vars.year)
    if testing == True:
        vars.solar_system.engine_settings(force_box, boxes, part_consumed_box, initial_fuel_mass, \
        t-t0*vars.year, planet_position_t0 + nt.rotate(vars.radius_AU[0]*unitvector, theta), t0)
        vars.solar_system.mass_needed_launch(final_position, test = False)

    launch_pos = planet_position_t0 + nt.rotate(vars.radius_AU[0]*unitvector, theta)
    #new_mass, dvv = boost(force, fuel_consumption, mass, satellite_mass, 1000000, 0.00001)
    #print('DELTA V', dvv*vars.year/vars.AU_tall)
    #print('MASS', new_mass)
    return t_AU, final_position, final_velocity, mass/vars.solmasse, fuel_mass/vars.solmasse, phi, launch_pos

def test(testing = False):
    t_load = np.load('saved/saved_orbits/launch_resolution/time_onlysun.npy')
    x_load = np.load('saved/saved_orbits/launch_resolution/pos_onlysun.npy')
    v_load = np.load('saved/saved_orbits/launch_resolution/vel_onlysun.npy')
    dt = t_load[1] - t_load[0]
    #print(planet_pos_t0)
    theta = 1/2*np.pi
    t0 = 0
    t, pos, vel, mass, fuel_mass, phi, launch_pos = launch(t_load, x_load, v_load, t0, theta, False)
    t, pos, vel, mass, fuel_mass, phi, launch_pos = launch(t_load, x_load, v_load, t0, theta-phi, True)
    t1_index = int((t)/dt)
    # print(t*vars.year)
    # print(pos*vars.AU_tall)
    # print(vars.radius[0]*1000)
    # print(vel*vars.AU_tall/vars.year)
    # print(mass*vars.solmasse)
    # print(fuel_mass*vars.solmasse)
    if testing == True:
        #index = min(range(len(t_load)), key=lambda i: abs(t_load[i]-t_AU))
        #print('INDEX', t1_index)
        measured_position = p4.position_from_objects(t1_index, vars.solar_system.analyse_distances(), x_load)
        print('position after launch from part4', measured_position)
        print('position from launch', pos)

        measured_velocity = p4.velocity_from_stars(vars.solar_system.measure_doppler_shifts())
        print('velocity after launch from part4', measured_velocity)
        print('velocity from launch.py', vel)
        print('position error', (measured_position-pos))
        print('velocity error', (measured_velocity-vel))
        print('Final Angle =', )
        vars.solar_system.take_picture()
        from PIL import Image
        find_orient = Image.open('find_orient.png')
        find_orient2 = np.array(find_orient)
        measured_angle = p4.find_angle(np.array(find_orient2))
        vars.solar_system.manual_orientation(measured_angle, measured_velocity, measured_position)


if __name__ == '__main__':
    #force_box = np.load('saved/engine/force_box.npy')
    #fuel_consumed_box_per_sec = np.load('saved/engine/fuel_consumed_box_per_sec.npy')
    #launch(force_box, fuel_consumed_box_per_sec, testing = True)
    test(True)
    pass

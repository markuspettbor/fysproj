import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from PIL import Image
import variables as vars
import numtools as nt
import orbit_tools as ot
import launch
import sys, os
import part6
import part7

# This code is ours; warning, this program creates 1000 files (orients in the current directory)

def find_launch_time(time, tol, x, v, mass, m_sat, target_indx, host_indx):
    '''
    Attempts to find launch time for satellite of mass m_sat
    time is time vector, tol is allowed deviation in AU between ideal orbit
    and target. x,v are orbits of all planets, target_indx, host_indx are
    indices of target and host planets.
    Assumes sun is index 0.
    '''
    steps = len(time)
    m_t = np.append(mass, m_sat)
    return ot.trajectory(m_t, x, v, host_indx, -1, target_indx, 0, time, tol)


def min_dist_res():
    # Minimum resolution distance
    r = radius[1]
    p = 1000
    f = 70
    print('Resolution distance:',r*p/f)

def check_orients(num_orients):
    # Read num_orients orients, spit out orbital parameters
    pos_sat = np.zeros((num_orients, 2))
    vel_sat = np.zeros((num_orients, 2))
    pos_host = np.zeros((num_orients, 2))
    pos_target = np.zeros((num_orients, 2))
    vel_target = np.zeros((num_orients, 2))
    time_orient = np.zeros((num_orients, 1))

    for i in range(num_orients):
        name = 'orient' + str(i) + '.npy'
        all_data = np.load(name)
        time_orient[i] = np.array(all_data[0])
        pos_planet = all_data[1]
        pos_vel_sat = all_data[2]
        pos_host[i] = pos_planet[:, host-1]
        pos_target[i] = pos_planet[:, target-1]
        pos_sat[i] = pos_vel_sat[:2]
        vel_sat[i] = pos_vel_sat[2:]
    return pos_sat, vel_sat, pos_host, pos_target, time_orient[:,0]

def get_dev(t_orient, optx, optv, x_sat, v_sat, nums, topt):
    # Find deviation from optimal orbit with optimal position optx, optimal velocity optv
    devx = np.zeros((nums-1))
    devv = np.copy(devx)
    dv = np.zeros((nums-1, 2))
    vec_between = np.copy(dv)
    for ii in range(nums-1):
        arg1 = np.argmin(np.abs(topt - t_orient[ii]))
        arg2 = np.argmin(np.abs(topt - t_orient[ii + 1]))
        devx[ii] =  nt.norm(x_sat[ii]-optx[arg1])
        devv[ii] =  nt.norm(v_sat[ii]-optv[arg1])
        dv[ii] = -v_sat[ii] + optv[arg1]
        vec_between[ii] = nt.unit_vector(-x_sat[ii] + optx[arg1])
    return devx, devv, dv, vec_between

def interpify(x1, t_orient):
    # Interpolate orbit
    xsx = interp1d(t_orient, x1[:,0])
    xsy = interp1d(t_orient, x1[:,1])
    x_sat_interp = lambda t: np.array([xsx(t), xsy(t)]).transpose()
    return x_sat_interp

def add_command(filename, t_of_boost, boost, command = 'boost', angle = [np.pi/2,np.pi/2], x = [0,0]):
    # Handle creation and changing command files for lander and satellite for MCast
    if command == 'createsat':
        with open(filename, 'w') as f:
            bo_str = 'launch\n'
    elif command == 'createlander':
        with open(filename, 'w') as f:
            bo_str = 'init\n'
    with open(filename, 'r') as f:
        lines = f.readlines()
    boost_ind = 1 # Launch at first pos
    if len(lines) > 1:
        for i in range(len(lines) - 1):
            next = lines[i+1].split()[1]
            if len(lines[i].split()) > 1:
                time = lines[i].split()[1]
                if float(time) <= t_of_boost and float(next) >= t_of_boost:
                    boost_ind = i + 1
        if t_of_boost > float(next):
            boost_ind = len(lines) + 1
    if command == 'boost':
        bo_str = 'boost '+str(t_of_boost)+' '+ str(boost[0])+' '+str(boost[1])+'\n'
    elif command == 'orient':
        bo_str = 'orient ' + str(t_of_boost) +'\n'
    elif command == 'picture':
        bo_str = 'picture ' + str(t_of_boost) +' '+ str(boost[0]) +' '+ str(boost[1])+' ' + '0' +' '+ '0' +' ' + ' 1' + '\n'
    elif command == 'video':
        bo_str = 'video ' + str(t_of_boost) + ' ' + str(angle[0]) + ' ' + str(angle[1]) + '\n'
    elif command == 'video_focus_on_planet':
        bo_str = 'video ' + str(t_of_boost) + ' ' + ' 1' + '\n'
    elif command == 'launchlander':
        bo_str = 'launchlander '+str(t_of_boost)+' '+ str(boost[0])+' '+str(boost[1])+' '+str(0)+'\n'
    elif command == 'boostlander':
        bo_str = 'boost '+str(t_of_boost)+' '+ str(boost[0])+' '+str(boost[1])+' '+str(0)+'\n'
    elif command == 'parachute':
        bo_str = 'parachute '+str(t_of_boost)+' '+ str(boost)+'\n'
    elif command == 'remove':
        bo_str = ''
        lines[boost_ind - 1] = ''
    lines.insert(boost_ind, bo_str)
    with open(filename, 'w') as f:
        for line in lines:
            print(line, file = f, end = '')

def interp_launch(filename):
    sys.stdout = open(os.devnull, 'w')
    solar_system.engine_settings(force_per_box, n_boxes, n_particles_sec_box,\
    initial_fuel, launch_dur, launch_pos, t_launch)
    solar_system.mass_needed_launch(fin_pos)
    solar_system.send_satellite(filename)
    sys.stdout = sys.__stdout__

m_star = vars.m_star
m_planets = vars.m
radius = vars.radius*1000/vars.AU_tall # Planet radii given in AUs
m_sat = vars.satellite/vars.solar_mass
x0 = np.array([[x0, y0] for x0, y0 in zip(vars.x0, vars.y0)])
v0 = np.array([[vx0, vy0] for vx0, vy0 in zip(vars.vx0, vars.vy0)])
x0 = np.concatenate((np.zeros((1,2)), x0))  # Set sun initial conditions
v0 = np.concatenate((np.zeros((1,2)), v0))
mass = np.append(m_star, m_planets) # Init solar system values

host = 1
target = 2 # Planet 1 in Mcast language
t0 = 0
t1 = 0.6
steps = 200000
time = np.linspace(t0, t1, steps)
#x, v = ot.patched_conic_orbits(time, mass, x0, v0)
#np.save('saved/saved_params/xx_sol2.npy', x)
#np.save('saved/saved_params/vv_sol2.npy', v)
x = np.load('saved/saved_params/xx_sol2.npy')
v = np.load('saved/saved_params/vv_sol2.npy')
tol = 5e-5
#t_launch_est, t_cept = find_launch_time(time, tol, x, v,\
#                                          mass, m_sat, target, host)
#np.save('saved/saved_params/t_launch1_est.npy', t_launch_est) # Launch window selector
#np.save('saved/saved_params/t_cept1.npy', t_cept)
t_launch_est = np.load('saved/saved_params/t_launch1_est.npy')
t_cept = np.load('saved/saved_params/t_cept1.npy')
t_launch = t_launch_est[0]
t_intercept = t_cept[0]
launch_indx = np.argmin(np.abs(t_launch - time))
intercept_indx = np.argmin((np.abs(t_intercept - time))) # Predicted index of intercept

print('Launch window selected at t=', t_launch)
print('Estimated time of intercept: t=', t_intercept)
x = x.transpose(); v = v.transpose()
fin_t, fin_pos, fin_vel, sat_mass, fuel_rem, angle, launch_pos = launch.launch(time, x, v, t_launch, testing = False)
x = x.transpose(); v = v.transpose()
force_per_box, n_boxes, n_particles_sec_box, initial_fuel, launch_dur = launch.get_engine_settings(t_launch_est[0], fin_t)
# Actually launch the thing
from ast2000solarsystem import AST2000SolarSystem
user = 'markusbpkjetilmg'
seed = AST2000SolarSystem.get_seed(user)
solar_system = AST2000SolarSystem(seed)
solar_system.engine_settings(force_per_box, n_boxes, n_particles_sec_box,\
initial_fuel, launch_dur, launch_pos, t_launch)
solar_system.mass_needed_launch(fin_pos)
launched_indx = np.argmin(np.abs(time-fin_t))

import part4 as p4
x = x.transpose(); v = v.transpose()
#Manual orientation
indx = launched_indx
measured_position = p4.position_from_objects(indx, solar_system.analyse_distances(), x)
measured_velocity = p4.velocity_from_stars(solar_system.measure_doppler_shifts())
solar_system.take_picture()
find_orient = Image.open('find_orient.png')
find_orient2 = np.array(find_orient)
measured_angle = p4.find_angle(np.array(find_orient2))
solar_system.manual_orientation(measured_angle, measured_velocity, measured_position)
x = x.transpose(); v = v.transpose()

topt = time[launch_indx]
t2 = time[launched_indx]
time2 = np.linspace(t2, t1, 150000) # New time vector, for optimal orbit from launch to end, high res

x0_sat = x[launch_indx, host]
v0_sat = v[launch_indx, host]
semimajor = (nt.norm(x0_sat) + nt.norm(x[intercept_indx, target]))/2

#x, v = ot.patched_conic_orbits(time2, mass, x[launched_indx], v[launched_indx])
#np.save('saved/saved_params/xx_sol4.npy', x)
#np.save('saved/saved_params/vv_sol4.npy', v)

x = np.load('saved/saved_params/xx_sol4.npy')
v = np.load('saved/saved_params/vv_sol4.npy')

dv = np.zeros((len(time2), 2))
launch_indx = 0
launched_indx = 0
intercept_indx = np.argmin(np.abs(t_intercept-time2))
dv_opt = ot.vis_viva(mass[0], nt.norm(x0_sat), semimajor)
deviation = np.pi/2 - nt.angle_between(v0_sat, x0_sat)
v0_opt = dv_opt*nt.rotate(nt.unit_vector(v0_sat), deviation)

acc = lambda r, t: ot.gravity(m_sat, m_star, r)/m_sat
t_optvec = np.linspace(topt, t2, 1000)
#a, b = ot.orbit(x0_sat, v0_opt, acc, t_optvec) # This only makes sure that the optimal orbit and the launch with mcast actually align.
#a = a.transpose(); b = b.transpose()
#opt_orb, opt_vel = ot.orbit(a[-1], b[-1], acc, time2) # Calculate optimal orbit
#np.save('saved/saved_params/xopt1.npy', opt_orb)
#np.save('saved/saved_params/vopt1.npy', opt_vel)
opt_orb = np.load('saved/saved_params/xopt1.npy')
opt_vel = np.load('saved/saved_params/vopt1.npy')
opt_orb = opt_orb.transpose()
opt_vel = opt_vel.transpose()


min_dist_res() # Minimum resolution distance
nums = 1000 # Number of orients. Warning: Will create 1000 files
t_start = time2[0]
t = np.linspace(t_start, time2[-1], nums)
print('time2[-1]', time2[-1])
add_command('satCommands.txt', 0, 0, command = 'createsat')
for tt in t:
    add_command('satCommands.txt', tt, 0, command = 'orient')
interp_launch('satCommands.txt') # Initial launch, only orients

x1, v1, p1, p2, t_orient = check_orients(nums) # Find orbital params
pint = interpify(p1, t_orient)
xs = interpify(x1, t_orient) # Interpolated orbit
vs = interpify(v1, t_orient)
req_boost_dist = ot.grav_influence(m_star, mass[host], xs(time2))[0] # Safe to boost
dist_to_host = nt.norm(xs(time2) - pint(time2), ax = 1)
first_boost = np.where(dist_to_host > 10*req_boost_dist)[0][0] # Even safer to boost
boost_time = time2[first_boost]
diff = np.zeros(2)

dt = t[1] - t[0]
min_altitude = radius[1]
max_altitude = 8e-5 # If the closest approach is less than this, stop optimizing
x_target = interpify(p2, t_orient)
num_boost = 10 # Number of corrections to make toward optimal orbit

dv_transfer = 0 # Delta v for transfer maneuver

for j in range(first_boost, first_boost + 1):
    # Orbital transfer optimization
    diff = np.zeros([num_boost, 2])
    boost_time_save = np.zeros(num_boost)
    for i in range(num_boost):
        boost_time = time2[j]
        vs = interpify(v1, t_orient)
        xs = interpify(x1, t_orient)
        dist_to_target = nt.norm(x_target(time2) - xs(time2), ax = 1)
        print('min(dist_to_target)', min(dist_to_target))
        if min(dist_to_target) > min_altitude and min(dist_to_target) < max_altitude:
            break
        diffv = opt_vel - vs(time2) # deviation in velocity
        diffx = opt_orb - xs(time2) # dev in position
        diff[i] = diffv[j] + diffx[j]*100 # weight position more
        print('DV', diffv[j])
        print('DX', diffx[j])
        dv_transfer += nt.norm(diff[i])
        add_command('satCommands2.txt', 0, 0, command = 'createsat')
        for tt in t_orient:
            add_command('satCommands2.txt', tt, 0, command = 'orient')
        boost_time_save[i] = boost_time
        for k in range(i+1):
            add_command('satCommands2.txt', boost_time_save[k], diff[k, :])
            num_boost_end = k
        interp_launch('satCommands2.txt') # launch, with more optimized orbit
        x1, v1, p1, p2, t_orient = check_orients(nums)
        j = int(j + first_boost*0.5)
        print('Time of optimization boost:', boost_time)

opt_transfer_boost = diff
xs = interpify(x1, t_orient) # Find sat orbit after optimization
vs = interpify(v1, t_orient)
dist_to_target = nt.norm(x_target(time2) - xs(time2), ax = 1) # Closest approach
p2 = interpify(p2, t_orient)
altitude = min(dist_to_target)
print('altitude', altitude)
print('radius[target-1]', radius[target -1])
periapsis = 1.5e-5 # If more than this value, we miss terribly, sadly. Value is inside planet
semi = (periapsis + altitude)/2 # Semimajor axis of elliptical orbit
print('periapsis', periapsis)
inject_point = np.argmin(dist_to_target) # When to inject
inject_time = time2[inject_point]
print('Time of injection:', inject_time)
v_target = np.gradient(p2(time2), axis = 0)/(time2[1]-time2[0]) # Approximate vel. of target planet
inject_vec = nt.unit_vector(x_target(inject_time) - xs(inject_time)) # Injection maneuver vector direction
orbital_vel = ot.vis_viva(mass[target], altitude, semi) # velocity at apoapsis
inject_vec = nt.rotate(inject_vec, -np.pi/2)*orbital_vel - vs(inject_time) + v_target[inject_point] # injection maneuver vector
kplanetgrav_at_inject = mass[target]/m_star*nt.norm(xs(inject_time))**2/dist_to_target[inject_point]**2
print('Planet gravity k = ', kplanetgrav_at_inject, ' stronger than star')
print('dv used, transfer:', dv_transfer)

nums = 1000 # 1000 orients
t_inject = np.linspace(time2[inject_point], time2[inject_point] + 0.005, nums)
dt = t_inject[1] - t_inject[0]
add_command('satCommands3.txt', 0, 0, command = 'createsat')
for tt in t_inject:
    add_command('satCommands3.txt', tt, 0, command = 'orient')
for ii in range(num_boost_end + 1):
    add_command('satCommands3.txt', boost_time_save[ii], opt_transfer_boost[ii, :]) # Transfer orbit
add_command('satCommands3.txt', inject_time, inject_vec) # Injection maneuver
interp_launch('satCommands3.txt') # Perform launch with injection, and optimized transfer orbit
x1, v1, p1, p2, t_orient = check_orients(nums)
num_interp = 10000 # Interpolate numinterp points near planet
allowed = 0
allowed_upper = int(num_interp/4)
dv_inject = nt.norm(inject_vec)
print('dv used, injection:', dv_inject)
dv_circ = 0

for i in range(4):
    # Circularization maneuver
    t_interp = np.linspace(t_inject[0], t_inject[-1], num_interp)
    xs = interpify(x1, t_inject)
    vs = interpify(v1, t_inject)
    p2 = interpify(p2, t_inject)
    v_target = np.gradient(p2(t_interp), axis = 0)/(t_interp[1]-t_interp[0])
    # Find lowest point in orbit, and only allow boosts after the first boost
    # To not boost at last point, upper limit is also set.
    circ_point = allowed + np.argmin(nt.norm(xs(t_interp[allowed:allowed_upper]) - p2(t_interp[allowed:allowed_upper]), ax = 1))
    circ_time = t_interp[circ_point]
    print('circ_point', circ_point)
    print('circ_time', circ_time)
    circ_radius = nt.norm(xs(circ_time) - p2(circ_time))
    print('circ_radius', circ_radius)
    vec_between = nt.unit_vector(xs(circ_time) - p2(circ_time))
    vec_between = nt.rotate(vec_between, np.pi/2) # Need perpendicular boost vector
    circ_vel = ot.vis_viva(mass[target], circ_radius, circ_radius)
    circularize_vec = -vs(circ_time) +  circ_vel*vec_between + v_target[circ_point]
    add_command('satCommands3.txt', circ_time, circularize_vec)
    interp_launch('satCommands3.txt') # Launch into proper circular orbit, with all maneuvers
    allowed = int(circ_point + 1)
    print('allowed', allowed)
    dv_circ += nt.norm(circularize_vec)
    allowed_upper = int(allowed_upper + num_interp/4)
    x1, v1, p1, p2, t_orient = check_orients(nums)

print('dv used, circ', dv_circ)
print('CLOSEST APPROACH:', np.min(nt.norm(x1 - p2, ax = 1)))
np.save('timefororbparam', t_orient)

def save_data():
    x1, v1, p1, p2, t_orient = check_orients(nums) #x1 = possat, v1 = velsat, p1 = posplan0, p2 = posplan1
    pos = x1-p2
    vel = v1
    angle = np.pi/3
    data = np.array([t_orient, pos, vel, p2, angle, -2])
    np.save('saved/saved_orbits/data_to_lander.npy', data)

def landing(nums):
    #The landing sequence, works with part 6 and 7
    x1, v1, p1, p2, t_orient = check_orients(nums) #x1 = possat, v1 = velsat, p1 = posplan0, p2 = posplan1
    boost = 0.8 #new velocity will be 80% of old velocity, in order to fall into planet
    pos_timezero = np.array([3266305.43741153, 1775670.55418141])
    angle_timezero = np.arctan2(pos_timezero[1], pos_timezero[0])
    if angle_timezero < 0:
        angle_timezero += 2*np.pi
    print('ANGLE ZERO = ', angle_timezero, angle_timezero*180/np.pi)

    pos_landing = np.array([260.606610, -3808732.23]) # x,y position when picture was taken at time = 6000, want to land there
    angle_landing = np.arctan2(pos_landing[1], pos_landing[0])
    if angle_landing < 0:
        angle_landing += 2*np.pi
    print('ANGLE WANTED = ', angle_landing, angle_landing*180/np.pi)
    angle_testing = part7.new_angle(angle_landing, 6137, 0)
    print('angle testing =', angle_testing, angle_testing*180/np.pi)
    print('angle testing - 5 =', angle_testing - angle_timezero, (angle_testing - angle_timezero)*180/np.pi)
    print('angle wanted, adjusted', angle_landing - angle_timezero, (angle_landing - angle_timezero)*180/np.pi)
    pos = x1 - p2
    vel_p2 = (p2[-1]-p2[-3])/(t_orient[-1]-t_orient[-3])
    vel = v1[-2] - vel_p2
    print('TIME DIFF END', t_orient[-1]-t_orient[-3])
    #vel_int = interpify(vel, t_orient)
    index_to_p7 = -3
    print('TIME IN AU TO SIM =', t_orient [-2])
    print('POS TO SIMULATION =', pos[index_to_p7,:]*vars.AU_tall)
    print('VEL TO SIMULATION =', vel*vars.AU_tall/vars.year)
    time_parachute, time_boost, boost_velocity, time_landed, angle_eject = \
            part7.optimise_landing(pos[index_to_p7,:]*vars.AU_tall, \
            vel*vars.AU_tall/vars.year, angle_landing, boost, plotting = True)
    angle_vector = np.arctan2(pos[1], pos[0])
    index = np.argmin(np.abs(angle_vector - angle_eject))
    offset_for_pictures = 0
    boost_lander = boost_velocity # Select vel for lander
    print('adding launchlander command')
    dt = t_orient[1] - t_orient[0]
    add_command('landerCommands3.txt', 0,0, command = 'createlander')
    add_command('landerCommands3.txt', time_boost + 1e-8 + offset_for_pictures, boost_lander, command = 'launchlander')
    area = 25 # Area of Parachute
    add_command('landerCommands3.txt', time_parachute + offset_for_pictures, area, command = 'parachute')
    #add_command('landerCommands3.txt', 3, 0, angle = np.array([0, 0]), command = 'video_focus_on_planet') # Create landing video
    #add_command('landerCommands3.txt', time_landed + 5000, 0, angle = np.array([0, 0]), command = 'video_focus_on_planet')
    add_command('landerCommands3.txt', time_landed + 4999, 0, command = 'orient') #need to give the porgram some time it can end, or it ends too early

    solar_system.land_on_planet(1, 'landerCommands3.txt')#, dt = 1) #LAND ON PLANETS

#these are to try to get values of grater resolution for initial contionions for the landing sequence
add_command('satCommands3.txt', t_inject[-1]+1e-8, 0, command = 'orient')
add_command('satCommands3.txt', t_inject[-1]+2e-8, 0, command = 'orient')
add_command('satCommands3.txt', t_inject[-1]+3e-8, 0, command = 'orient')
add_command('satCommands3.txt', t_inject[-1]+4e-8, 0, command = 'orient')
add_command('satCommands3.txt', t_inject[-1]+5e-8, 0, command = 'orient')
added = 5
interp_launch('satCommands3.txt') # Final launch, to send lander

print('saving')
save_data()
print('landing')
landing(nums + added)

def make_video():
    # Only used for making pretty video
    #add_command('satCommands3.txt', 0, 0, command = 'createsat')
    #for ii in range(num_boost_end + 1):
    #    print('ii', ii)
    #    add_command('satCommands3.txt', boost_time_save[ii], opt_transfer_boost[ii]) # Transfer orbit
    #add_command('satCommands3.txt', inject_time, inject_vec) # Injection maneuver
    diff_time = t_orient[-1] - inject_time
    add_command('satCommands3.txt', inject_time - diff_time/2 + 1e-8, 0, command = 'video_focus_on_planet') # movie start
    add_command('satCommands3.txt', t_orient[-1] - 1e-8, 0, command = 'video_focus_on_planet') # movie stop

    interp_launch('satCommands3.txt')
#make_video()

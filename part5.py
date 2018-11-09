import numpy as np
import matplotlib.pyplot as plt
import variables as vars
import numtools as nt
import orbit_tools as ot
import launch
import sys, os
import part6_kjetil as part6
import part7

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
    return ot.trajectory(m_t, x, v, host_indx, -1, target_indx, 0, time, True, tol)

m_star = vars.m_star
m_planets = vars.m
radius = vars.radius*1000/vars.AU_tall # Planet radii given in AUs
m_sat = vars.satellite/vars.solar_mass
x0 = np.array([[x0, y0] for x0, y0 in zip(vars.x0, vars.y0)])
v0 = np.array([[vx0, vy0] for vx0, vy0 in zip(vars.vx0, vars.vy0)])
x0 = np.concatenate((np.zeros((1,2)), x0))  # Set sun initial conditions
v0 = np.concatenate((np.zeros((1,2)), v0))
mass = np.append(m_star, m_planets)

host = 1
target = 2
t0 = 0
t1 = 0.6 # Heuristic value
steps = 200000
time = np.linspace(t0, t1, steps)
#x, v = ot.patched_conic_orbits(time, mass, x0, v0)
#np.save('saved/saved_params/xx_sol2.npy', x)
#np.save('saved/saved_params/vv_sol2.npy', v)
x = np.load('saved/saved_params/xx_sol2.npy')
v = np.load('saved/saved_params/vv_sol2.npy')
tol = 5e-5
#dv_est, t_launch_est, t_cept, semimajor = find_launch_time(time, tol, x, v,\
#                                           mass, m_sat, target, host)
#np.save('saved/saved_params/dv_est.npy', dv_est)
#np.save('saved/saved_params/t_launch_est.npy', t_launch_est)
#np.save('saved/saved_params/t_cept.npy', t_cept)
dv_est = np.load('saved/saved_params/dv_est.npy')
t_launch_est = np.load('saved/saved_params/t_launch_est.npy')
t_cept = np.load('saved/saved_params/t_cept.npy')

t_launch = t_launch_est[-1]
t_intercept = t_cept[-1]
launch_indx = np.argmin(np.abs(t_launch-time))
intercept_indx = np.argmin((np.abs(t_intercept - time)))

x = x.transpose(); v = v.transpose()
fin_t, fin_pos, fin_vel, sat_mass, fuel_rem, angle, launch_pos = launch.launch(time, x, v, t_launch, testing = False)
x = x.transpose(); v = v.transpose()
force_per_box, n_boxes, n_particles_sec_box, initial_fuel, launch_dur = launch.get_engine_settings(t_launch_est[-1], fin_t)
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
from PIL import Image
find_orient = Image.open('find_orient.png')
find_orient2 = np.array(find_orient)
measured_angle = p4.find_angle(np.array(find_orient2))
solar_system.manual_orientation(measured_angle, measured_velocity, measured_position)
x = x.transpose(); v = v.transpose()

t2 = time[launched_indx]

time2 = np.linspace(t2, t1, 150000)
'''
x, v = ot.patched_conic_orbits(time2, mass, x[launched_indx], v[launched_indx])
np.save('saved/saved_params/xx_sol3.npy', x)
np.save('saved/saved_params/vv_sol3.npy', v)
'''
x = np.load('saved/saved_params/xx_sol3.npy')
v = np.load('saved/saved_params/vv_sol3.npy')
launch_indx = 0
launched_indx = 0

dv = np.zeros((len(time2), 2))

x0_sat = fin_pos
v0_sat = fin_vel

intercept_indx = np.argmin(np.abs(t_intercept-time2))
semimajor = (nt.norm(x0_sat) + nt.norm(x[intercept_indx, target]))/2
dv_opt = ot.vis_viva(mass[0], nt.norm(x0_sat), semimajor)
deviation = np.pi/2 - nt.angle_between(v0_sat, x0_sat)
v0_opt = dv_opt*nt.rotate(nt.unit_vector(v0_sat), deviation)
'''
acc = lambda r, t: ot.gravity(m_sat, m_star, r)/m_sat
opt_orb, opt_vel = ot.orbit(x0_sat, v0_opt, acc, time2)
np.save('saved/saved_params/xopt.npy', opt_orb)
np.save('saved/saved_params/vopt.npy', opt_vel)
'''
opt_orb = np.load('saved/saved_params/xopt.npy')
opt_vel = np.load('saved/saved_params/vopt.npy')

opt_orb = opt_orb.transpose()
opt_vel = opt_vel.transpose()

#x_sat_init, v_sat_init, dv_used = ot.n_body_sat(x, mass, time2, dv, x0_sat, v0_sat, m_sat, opt_vel, opt_orb, 0, False)
'''
sphere_of_influence = ot.grav_influence(m_star, mass[host], x_sat_init)[0]
dist_to_host = nt.norm(x_sat_init - x[launched_indx:, host], ax = 1)
boost_time = np.argwhere(dist_to_host > 3*sphere_of_influence)[0,0]
closest = np.argmin(nt.norm(x_sat_init - x[:, target], ax = 1))
cdist = min(nt.norm(x_sat_init - x[:, target], ax = 1))
'''

def min_dist_res():
    r = radius[1]
    p = 1000
    f = 70
    print('Resolution distance:',r*p/f)
min_dist_res()

def real_launch(dv, tdv, numsteps, record = False):
    sys.stdout = open(os.devnull, 'w')
    solar_system.engine_settings(force_per_box, n_boxes, n_particles_sec_box,\
    initial_fuel, launch_dur, launch_pos, t_launch)
    solar_system.mass_needed_launch(fin_pos)
    with open('satCommands.txt', 'w') as f:
        print('launch', file = f)
        if record:
            print('video',str(tdv[0]), '1', file = f)
        for tt, dvv in zip(tdv, dv):
            print('orient', str(tt), file = f )
            if nt.norm(dvv) > 0:
                print('boost', str(tt), str(dvv[0]), str(dvv[1]), file = f)
        if record:
            print('video', str(tdv[-1]), '1', file = f)
            #print('video 0.55 1', file = f)
            #print('video 0.6 1', file = f)
    solar_system.send_satellite('satCommands.txt')
    sys.stdout = sys.__stdout__

def check_orients(num_orients):
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

from scipy.interpolate import interp1d

def interpify(x1, t_orient):
    xsx = interp1d(t_orient, x1[:,0])
    xsy = interp1d(t_orient, x1[:,1])
    x_sat_interp = lambda t: np.array([xsx(t), xsy(t)]).transpose()
    return x_sat_interp

nums = 200
t_start = time2[0]
t = np.linspace(t_start, time2[-1], nums)
diff = np.zeros((len(time2), 2))

real_launch(diff, t, nums)
x1, v1, p1, p2, t_orient = check_orients(nums)
pint = interpify(p1, t_orient)
p2int = interpify(p2, t_orient)
xs = interpify(x1, t_orient)
vs = interpify(v1, t_orient)

req_boost_dist = ot.grav_influence(m_star, mass[host], xs(time2))[0]
dist_to_host = nt.norm(xs(time2) - pint(time2), ax = 1)
first_boost = np.where(dist_to_host > 20*req_boost_dist)[0][0]
boost_time = time2[first_boost]
diff = np.zeros(2)


def interp_launch_commands(t, filename, launch = True, orient = True, record = False, r0 = 0, r1 = 1):
    with open(filename, 'w') as f:
        if launch:
            print('launch', file = f)
        if record:
            print('video', str(r0), '1', file = f)
        if orient:
            for ii in t:
                print('orient', str(ii), file = f)
        if record:
            print('video', str(r1), '1', file = f)

def add_command(filename, t_of_boost, boost, command = 'boost', angle = [np.pi/2,np.pi/2], x = [0,0]):
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
        bo_str = 'picture ' + str(t_of_boost) +' '+ str(angle[0]) +' '+ str(angle[1])+' ' + str(x[0]) +' '+ str(x[1]) +' ' + ' 1' + '\n'
    elif command == 'video':
        bo_str = 'video ' + str(t_of_boost) + ' ' + str(angle[0]) + ' ' + str(angle[1]) + '\n'
    elif command == 'video_focus_on_planet':
        bo_str = 'picture ' + str(t_of_boost) + ' ' + ' 1' + '\n'
    elif command == 'launchlander':
        bo_str = 'launchlander '+str(t_of_boost)+' '+ str(boost[0])+' '+str(boost[1])+' '+str(0)+'\n'
    elif command == 'boostlander':
        bo_str = 'boost '+str(t_of_boost)+' '+ str(boost[0])+' '+str(boost[1])+' '+str(0)+'\n'
    elif command == 'parachute':
        bo_str = 'parachute '+str(t_of_boost)+' '+ str(boost)+'\n'
    elif command == 'init':
        bo_str = 'init\n'
        boost_ind = 0
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

dt = t[1] - t[0]
min_altitude = radius[1]
max_altitude = 8e-5
x_target = interpify(p2, t_orient)

for j in range(first_boost, first_boost+1):
    for i in range(10):
        vs = interpify(v1, t_orient)
        xs = interpify(x1, t_orient)
        dist_to_target = nt.norm(x_target(time2) - xs(time2), ax = 1)
        print('min(dist_to_target)', min(dist_to_target))
        if min(dist_to_target) > min_altitude and min(dist_to_target) < max_altitude:
            break

        diffv = opt_vel - vs(time2)
        diffx = opt_orb - xs(time2)

        diff = diff + diffv[j] + diffx[j]
        interp_launch_commands(t_orient, 'satCommands2.txt')
        add_command('satCommands2.txt', boost_time, diff)
        interp_launch('satCommands2.txt')
        x1, v1, p1, p2, t_orient = check_orients(nums)

'''
plt.plot(x1[:,0], x1[:,1])
plt.scatter(x1[0,0], x1[0,1], c = 'm')
plt.axis('equal')
plt.show()
'''
opt_transfer_boost = diff
xs = interpify(x1, t_orient)
vs = interpify(v1, t_orient)
dist_to_target = nt.norm(x_target(time2) - xs(time2), ax = 1)
p2 = interpify(p2, t_orient)
altitude = min(dist_to_target)
print('altitude', altitude)
print('radius[target-1]', radius[target -1])
periapsis = 3.5e-5#radius[target - 1]
semi = (periapsis + altitude)/2
print('periapsis', periapsis)
inject_point = np.argmin(dist_to_target)
inject_time = time2[inject_point]
#v_target = v[inject_point, target]
v_target = np.gradient(p2(time2), axis = 0)/(time2[1]-time2[0])
inject_vec = nt.unit_vector(x_target(inject_time) - xs(inject_time))
orbital_vel = ot.vis_viva(mass[target], altitude, semi)#CHANGEDED HERE
inject_vec = nt.rotate(inject_vec, -np.pi/2)*orbital_vel*1 - vs(inject_time) + v_target[inject_point]

nums = 1000
t_inject = np.linspace(time2[inject_point], time2[inject_point] + 0.001, nums)
dt = t_inject[1] - t_inject[0]

interp_launch_commands(t_inject, 'satCommands3.txt', record = True, r0 = t_inject[0] - 0.001, r1 = t_inject[-1] + 0.01)
add_command('satCommands3.txt', boost_time, opt_transfer_boost) # Transfer orbit
add_command('satCommands3.txt', inject_time, inject_vec) # Injection maneuver
interp_launch('satCommands3.txt')
x1, v1, p1, p2, t_orient = check_orients(nums)

periapsis_indx = np.argmin(nt.norm(x1 - p2, ax = 1))

for i in range(1):
    t_interp = np.linspace(t_inject[0], t_inject[-1], 10000)
    xs = interpify(x1, t_inject)
    vs = interpify(v1, t_inject)
    p2 = interpify(p2, t_inject)

    circ_point = np.argmin(nt.norm(xs(t_interp) - p2(t_interp), ax = 1))
    circ_time = t_interp[circ_point]
    print('circ_point', circ_point)
    print('circ_time', circ_time)
    circ_radius = nt.norm(xs(circ_time) - p2(circ_time))
    print('circ_radius', circ_radius)
    vec_between = nt.unit_vector(xs(circ_time) - p2(circ_time))
    #angular_dev = np.pi/2 - nt.angle_between(vec_between, vs(circ_time))
    vec_between = nt.rotate(vec_between, np.pi/2)
    circ_vel = ot.vis_viva(mass[target], circ_radius, circ_radius)*1.021

    circularize_vec = -vs(circ_time) +  circ_vel*vec_between + v_target[np.argmin(np.abs(time2 - circ_time))]
    add_command('satCommands3.txt', circ_time, circularize_vec)
    interp_launch('satCommands3.txt')

    x1, v1, p1, p2, t_orient = check_orients(nums)

print('CLOSEST APPROACH:', np.min(nt.norm(x1 - p2, ax = 1)))


area = np.pi*radius[1]**2
print('area', area)
plt.scatter(0,0)
plt.plot(x1[:,0] - p2[:, 0], x1[:,1]- p2[:,1], '-k')
#plt.scatter(x1[0, 0], x1[0, 1], c = 'r')
plt.axis('equal')
#plt.show()



def plotting(nums):
    x1, v1, p1, p2, t_orient = check_orients(nums) #x1 = possat, v1 = velsat, p1 = posplan0, p2 = posplan1
    pi_vec = np.linspace(0, 2*np.pi, 1000)
    for theta in pi_vec:
        circle = part7.p2c_pos(np.array([radius[1]*1.27, theta]))
        plt.scatter(circle[0], circle[1], 0.1, 'k')
    pos = x1-p2
    vel = v1
    angle = np.pi/3 #TEMPORARY
    data = np.array([t_orient, pos, vel, p2, angle, 390])
    np.save('saved/saved_orbits/data_to_lander.npy', data)
    '''
    plt.scatter(pos[325,0], pos[325,1], c = 'r')
    plt.scatter(pos[350,0], pos[350,1], c = 'g')
    plt.scatter(pos[375,0], pos[375,1], c = 'b')
    plt.scatter(pos[400,0], pos[400,1], c = 'r')
    plt.scatter(pos[425,0], pos[425,1], c = 'g')
    plt.scatter(pos[450,0], pos[450,1], c = 'b')
    plt.scatter(pos[475,0], pos[475,1], c = 'r')
    '''
    plt.plot(pos[:,0], pos[:,1])
    plt.axis('equal')
    plt.show()

def landing(nums):
    x1, v1, p1, p2, t_orient = check_orients(nums) #x1 = possat, v1 = velsat, p1 = posplan0, p2 = posplan1
    boost = 0.8
    angle_landing = np.pi/3
    #x1_int = interpify(x1, t_orient)
    #p2_int = interpify(p2, t_orient)
    pos = x1 - p2
    #pos_int = x1_int - p2_int
    vel_p2 = np.gradient(p2, axis = 0)/(t_orient[1]-t_orient[0])
    vel = v1 - vel_p2
    #vel_int = interpify(vel, t_orient)
    #angle, time_parachute = part7.optimise_landing(t_orient, pos, vel, angle_landing, boost, 390)
    # TIME PARACHUTE ADDED CONVERTED TO YEARS AND ADDED TO TIME BOOST ?? ?
    angle_cheat = -0.26219225839919647
    angle_vector = np.arctan2(pos[1], pos[0])
    index = np.argmin(np.abs(angle_vector - angle_cheat))
    #index_int = np.argmin(np.abs(angle_vector - angle_cheat))
    boost_slow_time = t_orient[index]*vars.year
    boost_slow_velocity = -vel[index] + vel[index]*boost
    #COMMAND IS LAUNCHLANDER WITH TIME, VEL IS RELATIVE TO SATELITE (1-boost)*(-1*vel
    boost_lander = -vel[-1]*0.2*vars.AU_tall/vars.year
    print('adding launchlander command')
    dt = t_orient[1] - t_orient[0]
    interp_launch_commands(False, 'landerCommands3.txt', False, False)
    add_command('landerCommands3.txt', 0, 0, command = 'init')
    add_command('landerCommands3.txt', 1, boost_lander, command = 'launchlander')
    area = 6.2
    time_parachute = 2
    add_command('landerCommands3.txt', time_parachute, area, command = 'parachute')
    add_command('landerCommands3.txt', 10, 0, angle = np.array([1, '']), command = 'video')
    add_command('landerCommands3.txt', 12*60*60, 0, angle = np.array([1, '']), command = 'video')


    solar_system.land_on_planet(1, 'landerCommands3.txt') #LAND ON PLANETS
time = np.linspace(0.595,0.596, nums)
angle = np.linspace(0,2*np.pi, nums)
#for tid, vinkel in zip(time,  angle):
#    add_command('satCommands3.txt', tid, 0, command = 'orient')
interp_launch('satCommands3.txt')
print('landing')
landing(nums)
print('plotting')
plotting(nums)

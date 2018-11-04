import numpy as np
import matplotlib.pyplot as plt
import variables as vars
import numtools as nt
import orbit_tools as ot
import launch
import sys, os

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

def interp_launch(t, boost_time, diff, record = False):
    sys.stdout = open(os.devnull, 'w')
    solar_system.engine_settings(force_per_box, n_boxes, n_particles_sec_box,\
    initial_fuel, launch_dur, launch_pos, t_launch)
    solar_system.mass_needed_launch(fin_pos)
    with open('satCommands2.txt', 'w') as f:
        print('launch', file = f)
        if record:
            print('video', str(t[0]), '1', file = f)
        for ii in t:
            print('orient', str(ii), file = f)
            if ii <= boost_time and ii + dt >= boost_time:
                print('boost', str(boost_time), str(diff[j,0]), str(diff[j, 1]), file = f)
        if record:
            print('video', str(t[-1]), '1', file = f)
    solar_system.send_satellite('satCommands2.txt')
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
        print(min(dist_to_target))
        if min(dist_to_target) > min_altitude and min(dist_to_target) < max_altitude:
            break

        diffv = opt_vel - vs(time2)
        diffx = opt_orb - xs(time2)

        diff[j] = diff[j] + diffv[j] + diffx[j]
        #interp_launch(t, boost_time, diff)
        x1, v1, p1, p2, t_orient = check_orients(nums)

#interp_launch(t, boost_time, diff, record = True)  # Capture intercept on video
solar_system.engine_settings(force_per_box, n_boxes, n_particles_sec_box,\
initial_fuel, launch_dur, launch_pos, t_launch)
solar_system.mass_needed_launch(fin_pos)
solar_system.send_satellite('satCommands2.txt')
x1, v1, p1, p2, t_orient = check_orients(nums)
xs = interpify(x1, t_orient)
vs = interpify(v1, t_orient)
dist_to_target = nt.norm(x_target(time2) - xs(time2), ax = 1)

altitude = min(dist_to_target)
print(altitude)
inject_point = np.argmin(dist_to_target)
inject_time = time2[inject_point]
v_target = v[inject_point, target]
inject_vec = nt.unit_vector(x_target(inject_time)- xs(inject_time))
orbital_vel = ot.vis_viva(mass[target], altitude, altitude)
inject_vec = nt.rotate(inject_vec, -np.pi/2)*orbital_vel - vs(inject_time) + v_target

print(inject_vec, inject_time)

#interp_launch(t, boost_times, diff)
solar_system.engine_settings(force_per_box, n_boxes, n_particles_sec_box,\
initial_fuel, launch_dur, launch_pos, t_launch)
solar_system.mass_needed_launch(fin_pos)
solar_system.send_satellite('satCommands3.txt')
x1, v1, p1, p2, t_orient = check_orients(nums)

plt.plot(x1[:,0], x1[:,1])#, p1[:,0], p1[:,1], p2[:,0], p2[:,1])
plt.scatter(x1[:,0], x1[:,1])#, p1[:,0], p1[:,1], p2[:,0], p2[:,1])
#real_launch(ddd*0, t, nums)
#x1, v1, p1, p2, t_orient = check_orients(nums)
plt.scatter(x1[:,0], x1[:,1], c ='r')#, p1[:,0], p1[:,1], p2[:,0], p2[:,1])
plt.plot(opt_orb[:,0], opt_orb[:,1], '-r')
plt.axis('equal')
plt.show()

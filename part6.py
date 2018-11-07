def interp_launch_commands(t, filename, orient = True, record = False):
    with open(filename, 'w') as f:
        print('launch', file = f)
        if record:
            print('video', str(t[0]), '1', file = f)
        for ii in t:
            if orient:
                print('orient', str(ii), file = f)
        if record:
            print('video', str(t[-1]), '1', file = f)

def add_boost(filename, t_of_boost, boost, dt):
    with open(filename, 'r') as f:
        lines = f.readlines()
    boost_ind = len(lines)
    for i in range(len(lines)):
        if len(lines[i].split()) > 1:
            time = lines[i].split()[1]
            if float(time) < t_of_boost and float(time) + dt > t_of_boost:
                boost_ind = i
    bo_str = 'boost '+str(t_of_boost)+' '+ str(boost[0])+' '+str(boost[1])+'\n'
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
        print(min(dist_to_target))
        if min(dist_to_target) > min_altitude and min(dist_to_target) < max_altitude:
            break

        diffv = opt_vel - vs(time2)
        diffx = opt_orb - xs(time2)

        diff[j] = diff[j] + diffv[j] + diffx[j]

        interp_launch_commands(t_orient, 'satCommands2.txt')
        add_boost('satCommands2.txt', boost_time, diff[j], dt)
        interp_launch('satCommands2.txt')
        x1, v1, p1, p2, t_orient = check_orients(nums)

'''

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
        interp_launch(t, boost_time, diff)
        x1, v1, p1, p2, t_orient = check_orients(nums)
'''

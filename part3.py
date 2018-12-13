import numpy as np
import orbit_tools as ot
import numtools as nt
import variables as vars
import launch
import matplotlib.pyplot as plt

def solar_panel_area(watt, effic, F):
    return watt/(F*effic)*np.pi/2
radius_star = vars.radius_star*1000
dist_AU = np.sqrt(vars.x0**2 + vars.y0**2)
sigma = vars.sbc
T = vars.temp
radius = vars.radius_normal_unit #radius of planet
dist = dist_AU*vars.AU_tall
F = sigma*T**4 #flux per square meter of sun
surface_sun = 4*np.pi*(radius_star)**2
surface_shells = 4*np.pi*dist**2
L_tot = F*surface_sun #total W from the sun

F_rec = L_tot/surface_shells #W/m²
F_rec = F*surface_sun/surface_shells
solar_panel_A = solar_panel_area(40, 0.12, F_rec)
planet_temperature = T*(radius_star**2/dist**2/4)**(1/4)
planet_temperature = (F_rec/4/sigma)**(1/4)
#planet_temperature = (1/4*T**4*surface_sun/surface_shells)**(1/4)
Temps = np.array([390, 260])
distanser = np.sqrt(T**4/Temps**4*radius_star**2/4)
print(distanser/vars.AU_tall)
print(dist)

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

m_star = vars.m_star
m_planets = vars.m
a = vars.a
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
t1 = 0.6
period  = ot.kepler3(m_star, m_planets, a)[0]
orbits = t1/period
stepsperorbit = 10000
steps = stepsperorbit*orbits
print('Steps:', steps)
time = np.linspace(t0, t1, steps)
x, v = ot.patched_conic_orbits(time, mass, x0, v0)
tol = 5e-5
t_launch_est, t_cept = find_launch_time(time, tol, x, v, mass, m_sat, target, host)
t_launch = t_launch_est[0]
t_intercept = t_cept[0]
launch_indx = np.argmin(np.abs(t_launch - time))
intercept_indx = np.argmin((np.abs(t_intercept - time)))
print('Launch window selected at t =', t_launch)

t2 = time[launch_indx]
time2 = np.linspace(t2, t1, steps)
x0_sat = x[launch_indx, host]
v0_sat = v[launch_indx, host]
semimajor = (nt.norm(x0_sat) + nt.norm(x[intercept_indx, target]))/2
x, v = ot.patched_conic_orbits(time2, mass, x[launch_indx], v[launch_indx]) # Find planetary orbits
dv = np.zeros((len(time2), 2))
launch_indx = 0
intercept_indx = np.argmin(np.abs(t_intercept-time2))
dv_opt = ot.vis_viva(mass[0], nt.norm(x0_sat), semimajor)
deviation = np.pi/2 - nt.angle_between(v0_sat, x0_sat) # Angular deviation for Hohmann transfer
print('Angular deviation for optimal orbit', deviation)
v0_opt = dv_opt*nt.rotate(nt.unit_vector(v0_sat), deviation)
acc = lambda r, t: ot.gravity(m_sat, m_star, r)/m_sat
opt_orb, opt_vel = ot.orbit(x0_sat, v0_opt, acc, time2)
opt_orb = opt_orb.transpose(); opt_vel = opt_vel.transpose()
dvv = np.zeros((len(time2), 2))
v_escape = ot.vis_viva(m_planets[0], radius[0], 1e20) # Escape velocity
print('Escape velocity', v_escape)
dvv[0] = v_escape*nt.unit_vector(v0_sat) -v0_sat
x0_sat = x0_sat + nt.unit_vector(v0_sat)*radius[0]
x_sat, v_sat, dv_used = ot.n_body_sat(x, mass, time2, dvv, x0_sat, v0_sat, m_sat)
print('Closest approach:', min(nt.norm(x[:, 2] - x_sat, ax = 1)))
opt_dev = nt.norm(x_sat - opt_orb, ax = 1)
boost_thresh = ot.grav_influence(m_star, m_planets[0], x_sat, 1)
dist_to_host = nt.norm(x[:, 1] - x_sat, ax = 1)
print('Safe to boost for dist to host:', dist_to_host[0])
boost_time = np.where(dist_to_host > boost_thresh)[0][0]
print('Time of boost:', time2[boost_time])
x_sat, v_sat, dv_used = ot.n_body_sat(x, mass, time2, dvv, x0_sat, v0_sat, m_sat, opt_vel, opt_orb, 0, True)

closest = np.argmin(nt.norm(x[:, 2]- x_sat, ax = 1))
print('Closest approach:',min(nt.norm(x[:, 2] - x_sat, ax = 1)))
print('Time of closest approach:', time2[closest])
print('Predicted time of intercept:', t_intercept)
print('Total dv required:', np.sum(np.abs(dv_used)))
'''
plt.plot(x_sat[:, 0], x_sat[:,1], '--k', opt_orb[:, 0], opt_orb[:,1], '-.k',linewidth = 0.8)
plt.plot( x[:, 2, 0], x[:, 2,1 ], 'k', x[:, 1, 0], x[:, 1,1 ], 'b', linewidth = 0.8)
plt.xlabel('x [AU]', size = 12); plt.ylabel('y [AU]', size = 12)
plt.legend(['Satellite Orbit', 'Optimal Orbit', 'Target Planet Orbit', 'Home Planet Orbit'], frameon = False)
plt.scatter(x0_sat[0], x0_sat[1], c = 'k')
plt.scatter(x[intercept_indx, 2, 0], x[intercept_indx, 2, 1], c = 'k')
plt.scatter(x_sat[intercept_indx, 0], x_sat[intercept_indx, 1], c = 'k')
plt.scatter(x_sat[closest, 0], x_sat[closest, 1], c = 'k')
plt.scatter(x[closest, 2, 0], x[closest, 2, 1], c = 'k')
plt.axis('equal')
plt.show()
'''
if __name__ == '__main__':
    #print('At a distance %.6f AU from the sun, the flux equals %.3f W/m², and the lander needs %.3f m² of solar panels to function' %(dist_AU[0], F_rec, solar_panel_A))
    for i in [4,0,1,6,5,3,2]:
        print('\nPlanet number %i' %i)
        print('Solar panel area  ', solar_panel_A[i])
        print('Planet temperature', planet_temperature[i]-273.15)

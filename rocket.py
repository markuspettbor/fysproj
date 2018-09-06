import numtools as nt

def boost(self, thrust, consumption, initial_mass, delta_speed_desired, dt = 0.01):
    mass = initial_mass
    delta_speed = 0; a = 0
    while delta_speed < delta_speed_desired:
        delta_speed, mass = nt.euler_fuel_consumption(delta_speed, mass, thrust, consumption, dt)
    return initial_mass - mass

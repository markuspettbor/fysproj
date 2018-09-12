import variables as var
import numpy as np
import numtools as nt
import matplotlib.pyplot as plt

m1 = var.m
m2 = var.m_star
radius = var.radius
index = np.argmax(m1*m2/(radius**2))
rp = np.array([var.x0[index],var.y0[index]])
print(rp)


#create_radial_velocity()
vr = np.linspace(0,14,500)
vr = np.sin(vr)
vel = nt.create_radial_velocity(vr, 100, 3/7*np.pi)
plt.plot(vel)
plt.show()

def create_light_curve():

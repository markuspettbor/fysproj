import numpy as np
import matplotlib.pyplot as plt
import variables as vars
import part3

def density(radius): #h = altitude
    h = radius - vars.radius_normal_unit[1]
    try:
        h = float(h)
        val = True
    except TypeError:
        h = np.array(h)
        val = False
        rho = np.zeros(h.shape)

    rho0 = vars.rho0[1]
    T0 = part3.planet_temperature[1]
    m = (N2O + H2O)/2
    gamma = 1.4

    k = vars.k
    Tint = T0/2
    c = (rho0*k*T0/m)**(1-gamma)*T0**gamma
    M = vars.m_normal_unit[1]
    r0 = vars.radius_normal_unit[1]
    a = k*c**(1/gamma)/m
    b = vars.G_SI*M/a * (gamma - 1)/gamma

    rho1 = ( (rho0*a)**((gamma-1)) - \
        b*(1/r0 - 1/(r0 + h)) \
        )**(1/(gamma-1)) / a

    rho0_intersection = (c/(Tint)**gamma)**(1/(1-gamma)) * m/(k*Tint)
    if val == True:
        A = (rho0**(gamma-1) - rho0_intersection**(gamma-1))*a**(gamma-1)

        heigth = 1 / (1/r0-A/b) - r0

        #rho2 = np.exp(vars.G_SI*M*m/(vars.k*T0/2)*(1/(r0+h+heigth) - 1/(r0 + heigth)))*rho0_intersection
        #print(rho2)
        if h <= heigth:
            rho = rho1
        else:
            rho = np.exp(vars.G_SI*M*m/(vars.k*T0/2)*(1/(r0+h) - 1/(r0 + heigth)))*rho0_intersection
    else:
        rho_int = np.abs(rho1 - rho0_intersection)
        indx = np.nanargmin(rho_int)
        heigth = h[indx]
        rho2 = np.exp(vars.G_SI*M*m/(vars.k*T0/2)*(1/(r0+h+heigth) - 1/(r0 + heigth)))*rho0_intersection
        rho = rho1
        rho[indx:] = rho2[:(int(len(h) - indx))] #slices the arrays to make the final densityprofile
    return rho
#-------------------
oxy = 15.9994/vars.mol/1000
hyd = 1/vars.mol/1000
car = 12.0107/vars.mol/1000
nit = 14.00674/vars.mol/1000

O2  = 2*oxy
H2O = 2*hyd + oxy
CO2 = car + 2*oxy
CH4 = car + 4*hyd
CO  = car + oxy
N2O = 2*nit + oxy
#--------------------

if __name__ == '__main__':

    heigth = 300000
    res = 1/1000
    h = np.linspace(0, heigth, heigth*res+1) + vars.radius_normal_unit[1]
    for i in h:
        rho = density(i)
        plt.scatter(i,rho)
    #h = 10000

    #rho = density(h) #import denne fila, hent ut rho (tetthet)

    #print(rho)
    #plt.plot(h, rho, '-k')
    plt.show()

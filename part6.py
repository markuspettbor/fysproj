import numpy as np
from numba import jit
import variables as vars
import matplotlib.pyplot as plt
import part3_flux
#x_txt = np.loadtxt('saved/atmosphere/spectrum_seed09_600nm_3000nm.txt').transpose()
#print('loaded')
#np.save('saved/atmosphere/spectrum.npy', x_txt)
#sigma_noise_txt = np.loadtxt('saved/atmosphere/sigma_noise.txt').transpose()
#print('loaded')
#np.save('saved/atmosphere/sigma_noise.npy', sigma_noise_txt)



#print('loaded')
#print(measured_spectrum.shape)
'''
def density(h, rho0, T0, mH, gamma): #h = altitude
    h = np.array(h)
    rho = np.zeros(h.shape)
    konst = ((rho0*vars.k*T0)/(mu*mH)*T0**(gamma/(1-gamma)))**(1-gamma)
    M = vars.m_normal_unit
    rho = ( (rho0*vars.k*konst**(1/gamma) / (mu*mH))**((gamma-1)/gamma) - \
        mu*mH*vars.GSI*M/(k*konst**(1/gamma)) * (gamma - 1)/gamma*(1/r0 - 1/(r0 + h)) \
        )**(1/(gamma-1)) * mu*mH/(k*konst**(1/gamma))
    return rho
'''

#@jit(nopython = True)
def best_fit(a, b, c, f, noise, lambda_vector):
    best_sum = 1e100
    count = 0
    for z in c:
        for y in b:
            for x in a:
                g = (x)*np.exp(-(lambda_vector - y)**2/(2*z**2))
                g2 = 1-g
                funk = f - g2
                summ = np.sum(funk**2/noise)
                #print('One more down, many more to go')
                #print(lambda_vector)
                if summ < best_sum:
                    #print('NEW BEST')
                    best_sum = summ
                    best_x = x
                    best_y = y
                    best_z = z
                    count += 1
    #print('Count bests', count)
    return best_x, best_y, best_z

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
    T0 = part3_flux.planet_temperature[1]
    m = (N2O + H2O)/2 # same as mu*m_H
    gamma = 1.4

    k = vars.k
    Tint = T0/2
    c = (rho0*k*T0/m)**(1-gamma)*T0**gamma
    M = vars.m_normal_unit[1]
    r0 = vars.radius_normal_unit[1]
    a = k*c**(1/gamma)/m
    b = vars.G_SI*M/a * (gamma - 1)/gamma

    rho1 = ((rho0*a)**((gamma-1)/1) - b*(1/r0 - 1/(r0 + h)))**(1/(gamma-1)) / a
    #not 100% sure if the expression is correct..
    #calculates the density for adiabatic gas
    #runtime warning happens because the expression becomes negative, and this is taken to to the power of ~2.5
    plt.plot(rho1)
    plt.show()

    rho0_intersection = (c/(Tint)**gamma)**(1/(1-gamma)) * m/(k*Tint)
    if val == True: #If the input is a single height
        A = (rho0**(gamma-1) - rho0_intersection**(gamma-1))*a**(gamma-1)
        height = 1 / (1/r0-A/b) - r0

        #rho2 = np.exp(vars.G_SI*M*m/(vars.k*T0/2)*(1/(r0+h+height) - 1/(r0 + height)))*rho0_intersection
        #print(rho2)
        if h <= height:
            rho = rho1
        else:
            rho = np.exp(vars.G_SI*M*m/(vars.k*T0/2)*(1/(r0+h) - 1/(r0 + height)))*rho0_intersection
    else: #if the input is an array of heights
        rho_int = np.abs(rho1 - rho0_intersection)
        indx = np.nanargmin(rho_int)
        height = h[indx]
        rho2 = np.exp(vars.G_SI*M*m/(vars.k*T0/2)*(1/(r0+h+height) - 1/(r0 + height)))*rho0_intersection
        #calculated the density for isothermal gas
        print('Height at which adiabatic becomes isothermal is', height/1000, 'km')
        rho = rho1
        rho[indx:] = rho2[:(int(len(h) - indx))]
        #slices the arrays to make the final densityprofile, of both adia and iso
    return rho

#sets some global variables used in this part
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

def find_gasses():
    #siify = np.array([1, 1])
    sigma_noise = np.load('saved/atmosphere/sigma_noise.npy')
    #sigma_noise = (sigma_noise.transpose() * siify).transpose()
    measured_spectrum = np.load('saved/atmosphere/spectrum.npy')
    #measured_spectrum = (measured_spectrum.transpose() * siify).transpose()

    vmax = 11000
    vel = np.linspace(-vmax, vmax, 221)#, 2*vmax/1000*2 + 1)
    Tmin = 150
    Tmax = 450
    temp = np.linspace(Tmin, Tmax, 101)#, (Tmax-Tmin)/10*2 + 1)
    Fmin = 0.7
    F = np.linspace(Fmin, 1, 61)#, (1-Fmin)*10*4 + 1)

    lambda_0s = np.array([632, 690, 760, 720, 820, 940, 1400, 1600, 1660, 2200, 2340, 2870])
    masses =             [O2,  O2,  O2,  H2O, H2O, H2O, CO2,  CO2,  CH4,  CH4,  CO,   N2O]

    best_flux = np.zeros(len(lambda_0s))
    best_lambda = np.zeros(len(lambda_0s))
    best_sigma = np.zeros(len(lambda_0s))

    fluxes = 1 - F
    gaussiums = np.zeros([len(lambda_0s), len(measured_spectrum[0])])

    for i in range(len(lambda_0s)):
        lam0 = lambda_0s[i]
        mass = masses[i]
        lambdas = lam0 + vel/vars.c*lam0
        sigmas =  np.sqrt(vars.k*temp/mass)*lam0/vars.c

        span = 5*max(sigmas) + max(lambdas) - min(lambdas)
        idx_low = (np.abs(measured_spectrum[0] - (lam0 - span/2))).argmin()
        idx_high = (np.abs(measured_spectrum[0] - (lam0 + span/2))).argmin()
        spectrum = measured_spectrum[0, idx_low:idx_high]
        spectrum_values = measured_spectrum[1, idx_low:idx_high]
        noise_vect = sigma_noise[1, idx_low:idx_high]
        best_flux[i], best_lambda[i], best_sigma[i] = best_fit(fluxes, lambdas, sigmas, spectrum_values, noise_vect, spectrum)
        print('Number %i of %i complete' %((i+1), int(len(lambda_0s))))
        plt.figure()
        plt.plot(spectrum, spectrum_values)
        plt.plot(spectrum, 1 - (best_flux[i])*np.exp(-(spectrum-best_lambda[i])**2/(2*best_sigma[i]**2)))
        plt.title('Lambda0 = %f' %(lam0))
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Flux')

        #plt.plot(gaussiums[i])
        #plt.show()
    print('best', best_lambda)
    print('lam0', lambda_0s)
    g_vel = (best_lambda - lambda_0s)*vars.c/lambda_0s
    g_temp = best_sigma**2*masses/vars.k*vars.c**2/lambda_0s**2
    g_flux = 1 - best_flux

    for Vel, Temp, Flux, Lam0 in zip(g_vel, g_temp, g_flux, lambda_0s):
        if Flux != 1:
            print('Lambda = ', Lam0)
            print('Vel = %f, Temp = %f, Flux = %f' % (Vel, Temp, Flux))

    plt.show()

def test_landing():
    height = 160000
    res = 1/10
    radius = np.linspace(0, height, height*res+1) + vars.radius_normal_unit[1]

    #for i in h:
    #    rho = density(i)
    #    plt.scatter(i,rho)

    rho = density(radius) #import denne fila, hent ut rho (tetthet)

    plt.plot(radius/1000 - vars.radius_normal_unit[1]/1000, rho, '-k')
    plt.xlabel('Height [km]')
    plt.ylabel('Atmospheric Density [atm]')
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    #find_gasses()
    test_landing()

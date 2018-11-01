import numpy as np
from numba import jit
import variables as vars
import matplotlib.pyplot as plt
#x_txt = np.loadtxt('saved/atmosphere/spectrum_seed09_600nm_3000nm.txt').transpose()
#print('loaded')
#np.save('saved/atmosphere/spectrum.npy', x_txt)
#sigma_noise_txt = np.loadtxt('saved/atmosphere/sigma_noise.txt').transpose()
#print('loaded')
#np.save('saved/atmosphere/sigma_noise.npy', sigma_noise_txt)

siify = np.array([1, 1])

sigma_noise = np.load('saved/atmosphere/sigma_noise.npy')
sigma_noise = (sigma_noise.transpose() * siify).transpose()
#print(sigma_noise.shape)
measured_spectrum = np.load('saved/atmosphere/spectrum.npy')

measured_spectrum = (measured_spectrum.transpose() * siify).transpose()

#print('loaded')
#print(measured_spectrum.shape)

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
vmax = 13000
vel = np.linspace(-vmax, vmax, 27)#, 2*vmax/1000*2 + 1)
Tmin = 150
Tmax = 450
temp = np.linspace(Tmin, Tmax, 11)#, (Tmax-Tmin)/10*2 + 1)
Fmin = 0.7
F = np.linspace(Fmin, 1, 26)#, (1-Fmin)*10*4 + 1)

lambda_0s = np.array([632, 690, 760, 720, 820, 940, 1400, 1600, 1660, 2200, 2340, 2870])
masses =    [O2,  O2,  O2,  H2O, H2O, H2O, CO2,  CO2,  CH4,  CH4,  CO,   N2O]

best_flux = np.zeros(len(lambda_0s))
best_lambda = np.zeros(len(lambda_0s))
best_sigma = np.zeros(len(lambda_0s))

fluxes = 1 - F
gaussiums = np.zeros([len(lambda_0s), len(measured_spectrum[0])])

print('3, 2, 1, GO!')

for i in range(len(lambda_0s)):
    lam0 = lambda_0s[i]
    mass = masses[i]
    lambdas = lam0 + vel/vars.c*lam0
    sigmas =  lam0*vars.k*temp/vars.c/mass

    span = 5*max(sigmas) + max(lambdas) - min(lambdas)
    idx_low = (np.abs(measured_spectrum[0] - (lam0 - span/2))).argmin()
    idx_high = (np.abs(measured_spectrum[0] - (lam0 + span/2))).argmin()
    spectrum = measured_spectrum[0, idx_low:idx_high]
    spectrum_values = measured_spectrum[1, idx_low:idx_high]
    noise_vect = sigma_noise[1, idx_low:idx_high]
    #print('span =', span)
    #print('shape =', spectrum.shape)

    #g = lambda flu, lam, sig: (1 - flu)*np.exp(-(measured_spectrum[0]-lam)**2/(2*sig**2))
    best_flux[i], best_lambda[i], best_sigma[i] = best_fit(fluxes, lambdas, sigmas, spectrum_values, noise_vect, spectrum)
    print('Number %i of %i complete' %((i+1), int(len(lambda_0s))))
    print('Flux = %f, Lambda = %f, Sigma = %f' %(best_flux[i], best_lambda[i], best_sigma[i]))
    #gaussiums[i] = (best_flux[i])*np.exp(-(spectrum-best_lambda[i])**2/(2*best_sigma[i]**2))
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
g_temp = best_sigma*vars.c*masses/vars.k/lambda_0s
g_flux = 1 - best_flux

for Vel, Temp, Flux, Lam0 in zip(g_vel, g_temp, g_flux, lambda_0s):
    if Flux != 1:
        print('Lambda = ', Lam0)
        print('Vel = %f, Temp = %f, Flux = %f' % (Vel, Temp, Flux))

plt.show()
'''
print('flux',best_flux)
print('lambda',best_lambda)
print('sigma',best_sigma)
print(measured_spectrum.shape)


print('GOOOOOOAL')

continium = np.ones(len(measured_spectrum[0]))
print(gaussiums.shape, 'SHAPE')
estimatium = continium - np.sum(gaussiums, axis = 0)
print(estimatium.shape, 'SHAPE')

velocity = (best_lambda - lambda_0s)*vars.c/lambda_0s
temperature = best_sigma*vars.c*masses/lambda_0s/vars.k
fluxes_final = best_flux

print('velocity', velocity)
print('temperature', temperature)
print('fluxes', fluxes_final)

plt.plot(measured_spectrum[0]*1e9, measured_spectrum[1], '-k')
plt.plot(measured_spectrum[0]*1e9, estimatium, '-r')
plt.show()

print('Thank you for playing!')
'''

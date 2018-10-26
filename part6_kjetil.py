import numpy as np
import variables as vars
import matplotlib.pyplot as plt
#x_txt = np.loadtxt('saved/atmosphere/spectrum_seed09_600nm_3000nm.txt').transpose()
#print('loaded')
#np.save('saved/atmosphere/spectrum.npy', x_txt)
#sigma_noise_txt = np.loadtxt('saved/atmosphere/sigma_noise.txt').transpose()
#print('loaded')
#np.save('saved/atmosphere/sigma_noise.npy', sigma_noise_txt)
sigma_noise = np.load('saved/atmosphere/sigma_noise.npy')
#print(sigma_noise.shape)
measured_spectrum = np.load('saved/atmosphere/spectrum.npy')
#print('loaded')
#print(measured_spectrum.shape)

#@jit(nopython = True) #HOW TO RUN JIT??
def best_fit(a, b, c, f, g, noise):
    '''x = ..y = ..z =...'''
    best_sum = 1e100
    count = 0
    for x in a:
        for y in b:
            for z in c:
                summ = np.sum(((f - g(x,y,z))/noise)**2)
                if summ < best_sum:
                    best_sum = summ
                    best_x = x
                    best_y = y
                    best_z = z
                    count += 1
    return best_x, best_y, best_z
    print(count)

oxy = 15.9994/vars.mol
hyd = 1/vars.mol
car = 12.0107/vars.mol
nit = 14.00674/vars.mol
O2  = 2*oxy
H2O = 2*hyd + oxy
CO2 = car + 2*oxy
CH4 = car + 4*hyd
CO  = car + oxy
N2O = 2*nit + oxy

vmax = 13000
vel = np.linspace(-vmax, vmax, 2*vmax/1000 + 1)
Tmin = 150
Tmax = 450
temp = np.linspace(Tmin, Tmax, (Tmax-Tmin)/10 + 1)
Fmin = 0.7
F = np.linspace(Fmin, 1, (1-Fmin)*10 + 1)

lambda_0s = np.array([632, 690, 760, 720, 820, 940, 1400, 1600, 1660, 2200, 2340, 2870])*1e-9
masses =    [O2,  O2,  O2,  H2O, H2O, H2O, CO2,  CO2,  CH4,  CH4,  CO,   N2O]

best_flux = np.zeros(len(lambda_0s))
best_lambda = np.zeros(len(lambda_0s))
best_sigma = np.zeros(len(lambda_0s))

fluxes = 1 - F
gaussiums = np.zeros([len(lambda_0s), len(measured_spectrum[0])])

for i in range(len(lambda_0s)):
    lam0 = lambda_0s[i]
    mass = masses[i]
    lambdas = lam0 + vel/vars.c*lam0
    sigmas =  lam0/vars.k*temp/vars.c/mass
    g = lambda flu, lam, sig: (1 - flu)*np.exp(-(measured_spectrum[0]-lam)**2/(2*sig**2))
    best_flux[i], best_lambda[i], best_sigma[i] = best_fit(fluxes, lambdas, sigmas, measured_spectrum[1], g, sigma_noise[1])
    print('Number %i of %i complete' %((i+1), int(len(lambda_0s))))
    gaussiums[i] = g(best_flux[i], best_lambda[i], best_sigma[i])
    print(gaussiums[i].shape, 'SHAPE GAUSS')
    break
print('flux',best_flux)
print('lambda',best_lambda)
print('sigma',best_sigma)

continium = np.ones(len(measured_spectrum[0]))
print(gaussiums.shape, 'SHAPE')
estimatium = continium - np.sum(gaussiums)
print(estimatium.shape, 'SHAPE')
plt.plot(measured_spectrum[0], measured_spectrum[1], '-k')
plt.plot(measured_spectrum[0], estimatium, '-r')
plt.show()

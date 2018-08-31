import numpy as np
import prob_dist as pd
import matplotlib.pyplot as plt
import numtools as nt
'''
maxwell boltzmann
p(v) = (m/2pikT)**3/2*e**(-1/2*mv**2/kT)
--> black box --> sig = sqrt(kT/m), mu = 0

k =  1.380 648 52 x 10-23 J K-1
https://physics.nist.gov/cgi-bin/cuu/Value?k|search_for=boltzmann

Na = 6.022 140 857 x 1023 mol-1
https://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=avogadro


mass H_2 = 2.016 g/mol
https://pubchem.ncbi.nlm.nih.gov/compound/Hydrogen
'''
np.random.seed(69)
'''
molmassH2 = 2.016 #g/mol
mol = 6.022140857e23 #1/mol
k = 1.38064852e-23
T = 1e4
m = molmassH2 / mol / 1000
N = 10**5
sig = np.sqrt(k*T/m)
mu = 0
x = np.random.normal(mu, sig, (1,N))[0]
y = np.random.normal(mu, sig, (1,N))
z = np.random.normal(mu, sig, (1,N))

x = np.linspace(-2.5, 2.5, 10000)*10**4
gauss = pd.NormDist(mu, sig)
plt.plot(x, gauss.density_func(x))
plt.title('linspaces')
plt.show()
'''
# Alternative implementation:
x = np.linspace(-2.5, 2.5, 10000)*10**4
boltzmax = pd.BoltzmannMaxwell(T = 1e4, N = 10e5)
plt.plot(x, boltzmax(x))
plt.show()



# plt.hist(x, bins = 690, range = (-25000, 25000), color = 'k')
# plt.xlabel('Velocity in x-Direction')
# plt.ylabel('Number of Particles')
# plt.show()

# x_y, x_x = np.histogram(x,)
# x_where = np.where(np.logical_and(x_x > -30000, x_x < 30000))[0]
# print(x_where)
# x_integrate = x_y[x_where[0] : x_where[-1]]
# print(x_integrate)
# integ = np.trapz(x_integrate, x_x[x_where[0] : x_where[-1]])
# print(integ)
# print(x_where[-1])
# print(x_integrate)

# x_sort = np.sort(x)[0]
# x_where = np.where(np.logical_and(x_sort >= 5000, x_sort <= 30000))[0]
# print(x_sort)
# #funk = x_sort
# #x = x_where
# plt.plot(x_sort)
# plt.show()
# integral = nt.integrate(lambda x: x_sort[x_where[0]:x_where[-1]], x_where)
# print(integral)

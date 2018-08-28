import numpy as np
import prob_dist as pd
import matplotlib.pyplot as plt

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
molmassH2 = 2.016 #g/mol
mol = 6.022140857e23 #1/mol

k = 1.38064852e-23
T = 10*4
m = molmassH2 / mol / 1000
N = 10**5
sig = (k*T/m)**0.5
mu = 0
x = np.random.normal(mu, sig, (N,1))
y = np.random.normal(mu, sig, (N,1))
z = np.random.normal(mu, sig, (N,1))

plt.hist(x, bins = 69, range = (-2500, 2500), color = 'k')
plt.xlabel('Velocity in x-Direction')
plt.ylabel('Number of Particles')
plt.show()

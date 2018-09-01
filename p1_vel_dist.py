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
x = np.linspace(-2.5, 2.5, 10000)*10**4
boltzmax = pd.BoltzmannMaxwell(T = 1e4, N = 10e5)
plt.plot(x, boltzmax(x))
plt.show()

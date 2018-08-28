import numpy as np
import prob_dist as pd

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

k = 1.38e-23
T = 10*4
m =
N = 10**5
sig = (k*T/m)**0.5
mu = 0
np.random.normal(mu, sig, (1,N))

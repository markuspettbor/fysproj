import numpy as np
import prob_dist as pd


'''
maxwell boltzmann
p(v) = (m/2pikT)**3/2*e**(-1/2*mv**2/kT)
--> black box --> sig = sqrt(kT/m), mu = 0
'''

k = 1.38e-23
T = 10*4
m =
N = 10**5
sig = (k*T/m)**0.5
mu = 0
np.random.normal(mu, sig, (1,N))

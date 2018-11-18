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
import matplotlib

# Task 1
np.random.seed(10)
n = 10**5
x0 = np.linspace(-2.5, 2.5, 10000)*10**4
boltzmax = pd.BoltzmannMaxwell(T = 1e4, N = n)
plt.plot(x0, boltzmax(x0), '-k', linewidth = 0.8)
plt.xlabel('$v_x$ [m/s]', size = 12)
plt.ylabel('$P(v_x)$', size  = 12)
plt.title('Velocity Distribution, $P(v_x)$')
plt.show()

# Task 2
x1 = np.linspace(5, 30, 10000)*10**3
area = nt.integrate(lambda x: boltzmax(x1), x1, teacher_is_strict = False)
print('Integral of P(vx):', area, 'N*area:', n*area)

# Task 3
x2 = np.linspace(0, 3, 10000)*10**4
absolute_boltzmax = pd.AbsoluteBoltzmannMaxwell(T = 1e4, N = n)
plt.plot(x2, absolute_boltzmax.absolute_density(x2), '-k', linewidth = 0.8)
plt.xlabel('$v$ [m/s]', size = 12)
plt.ylabel('$P(v)$', size = 12)
plt.title('Absolute Velocity Distribution, $P(v)$')
plt.show()

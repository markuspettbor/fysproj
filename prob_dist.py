import numpy as np
import numtools as nt
import matplotlib.pyplot as plt

class ProbDist:
    def __init__(self, mu, sig):
        self.mu = mu
        self.sig = sig

    def density_func(self, x):
        # Density function of choice
        pass

    def getProb(self, x):
        return nt.integrate(self.density_func, x)

class NormDist(ProbDist):
    def density_func(self, x):
        return 1/(np.sqrt(2*np.pi)*self.sig)\
                *np.exp(-0.5*((x - self.mu)/self.sig)**2)

class BoltzmannMaxwell(NormDist):
    def __init__(self, T, N):
        mmH2 = 2.016 #g/mol
        mol = 6.022140857e23 #1/mol
        self.k = 1.38064852e-23 # Boltzmann constant
        self.m = mmH2/mol/1000
        self.T = T
        self.N = N
        self.sig = np.sqrt(self.k*self.T/self.m)
        self.mu = 0

    def __call__(self, v):
        return self.density_func(v)

    def distribution(self, size):
        return np.random.normal(self.mu, self.sig, size = size)

class AbsoluteBoltzmannMaxwell(BoltzmannMaxwell):
    def absolute_density(self, v):
        return 4*v**2*np.pi/(np.sqrt(2*np.pi)*self.sig)**2*self(v)

def test():
    mu = 100
    sigma = 5
    t = NormDist(mu, sigma)
    for i in range(1, 4):
        range_x = np.linspace(mu - i*sigma, mu + i*sigma, 1000)
        probability = t.getProb(range_x)
        print('Range [%.2f, %.2f] (%d sigma), with estimated probability %f.'
             %(min(range_x), max(range_x), i, probability))
    plt.plot(range_x, t.density_func(range_x))
    plt.show()

def test_bm():
    bolt = AbsoluteBoltzmannMaxwell(1000, 1000)
    v = np.linspace(0, 10, 1000)*10**3
    plt.plot(v, bolt.absolute_density(v))
    plt.show()


if __name__ == '__main__':
    test()
    test_bm()

import numpy as np
import numtools as nt

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

def test():
    mu = 100
    sigma = 5
    t = NormDist(mu, sigma)
    for i in range(1, 4):
        range_x = np.linspace(mu - i*sigma, mu + i*sigma, 1000)
        probability = t.getProb(range_x)
        print('Range [%.2f, %.2f] (%d sigma), with estimated probability %f.'
             %(min(range_x), max(range_x), i, probability))

if __name__ == '__main__':
    test()

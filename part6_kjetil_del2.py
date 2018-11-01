import numpy as np

T0_C = 5.054794035754696 #cel
T0 = T0_C + 273.15
P0 = 6.4304675739249424
gamma = 1.4
konst = P0**(1-gamma)*T0**gamma
print(konst)

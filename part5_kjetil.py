import numpy as np
stop = 3
x = np.linspace(0,stop,stop+1, dtype = 'int')
for i in x:
    for j in x[x != i]:
        for k in x[(x != i) * (x != j)]:
            print(i,j,k)

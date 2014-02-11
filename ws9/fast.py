#!/usr/bin/env python

import numpy as np
from timeit import timeit

times = []

# Function to read matrix data, check for invertibility of A and that the 
# sizes match, and return A and b.
def solve(i):
    j = i+1
    A = np.genfromtxt('LSE' + str(j) + '_m.dat')
    b = np.genfromtxt('LSE' + str(j) + '_bvec.dat')
    (a1, a2) = np.shape(A)
    b1 = len(b)
    if a1 != a2 or a1 != b1:
        print 'System cannot be solved. Dimension error.'
        return 0
    else:
        det = np.linalg.det(A)
        if det == 0:
            print 'System cannot be solved; A' + str(j) + ' not invertible.'
            return 0
        else:
            return (A,b)
            #return np.linalg.solve(A,b)

# I have since learned that linalg.solve already checks for singularity and
# non-squareness


for i in range(5):
    times.append(timeit("np.linalg.solve(A,b)", number=1, \
                 setup="import numpy as np; from fast import solve; (A,b)=solve(%d)" % i))

print times

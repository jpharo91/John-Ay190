#!/usr/bin/env python

import numpy as np

for i in range(5):
    j = i+1
    A = np.genfromtxt('LSE' + str(j) + '_m.dat')
    b = np.genfromtxt('LSE' + str(j) + '_bvec.dat')
    (a1, a2) = np.shape(A)
    b1 = len(b)
    print 'A' + str(j) + ': ' + str(np.shape(A))
    print 'b' + str(j) + ': ' + str(np.shape(b))
    if a1 != a2 or a1 != b1:
        print 'System ' + str(j) + ' cannot be solved. Dimension error.'
    else:
        det = np.linalg.det(A)
        if det == 0:
            print 'System cannot be solved; A' + str(j) + ' not invertible.'
        else:
            print 'System ' + str(j) + ' can be solved.'

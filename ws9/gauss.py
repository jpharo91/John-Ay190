#!/usr/bin/env python

import numpy as np
from timeit import timeit

# Function for Gauss elimination found from the glorious internet, specifically
# http://ine.scripts.mit.edu/blog/2011/05/gaussian-elimination-in-python/
def myGauss(m):
    #eliminate columns
    for col in range(len(m[0])):
        for row in range(col+1, len(m)):
            r = [(rowValue * (-(m[row][col] / m[col][col]))) for rowValue in m[col]]
            m[row] = [sum(pair) for pair in zip(m[row], r)]
    #now backsolve by substitution
    ans = []
    m.reverse() #makes it easier to backsolve
    for sol in range(len(m)):
            if sol == 0:
                ans.append(m[sol][-1] / m[sol][-2])
            else:
                inner = 0
                #substitute in all known coefficients
                for x in range(sol):
                    inner += (ans[x]*m[sol][-2-x])
                #the equation is now reduced to ax + b = c form
                #solve with (c - b) / a
                ans.append((m[sol][-1]-inner)/m[sol][-sol-2])
    ans.reverse()
    return ans


testMatrix = [[3.0, 2.0, -1.0, 1.0],[2.0, -2.0, 4.0, -2.0],[-1.0, 0.5, -1.0, 0.0]]

answer = myGauss(testMatrix)
for i in range(len(answer)):
    answer[i] = round(answer[i], 1)
print answer


answer2 = myGauss([[2.0,1.0,-1.0,8.0],
               [-3.0,-1.0,2.0,-11.0],
               [-2.0,1.0,2.0,-3.0]])
print answer2

times = []

# Function to read the matrix data and return an augmented matrix for use in
# Gaussian Elimination.
def data(i):
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
                b = [b]
                augment = np.concatenate((np.array(A),np.array(b).T), axis=1)
                return list(augment)

for i in range(5):
    times.append(timeit("myGauss(m)", number=1, \
             setup="from gauss import myGauss; from gauss import data; m=data(%d)" % i))
    print times

print times

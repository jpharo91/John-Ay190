#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl

# get and slice data
data = np.loadtxt('cepheids.txt', comments= '#')
t = data[:,0]
m = data[:,1]

# helper function to calculate fractions in L_nj
def l(x, j, k):
    return (x - k) / (j - k)

# function representing L_nj(x)
def L(x, j):
    product = 1
    k = [i for i in t if i != t[j]]  #list comprehension to prevent /0 error
    for q in k:
        product *= l(x, t[j], q)
    return product

# function to generate a Lagrange polynomial of degree n 
def p(n, x):
    total = 0
    for j in range (n+1):
        total += L(x, j)*m[j]
    return total

# generate x values and y values from interpolation
x = np.arange(0.0, 1.0, 0.001)
y = p(8, x)

# plot data and interpolation
p1, = pl.plot(t, m, 'ro', linewidth=2)
p2, = pl.plot(x, y, "b", linewidth=2)

# set x and y ranges
pl.xlim(min(t),max(t))
pl.ylim(min(y-.25),max(y*1.05))

# label the axes
pl.xlabel("Time")
pl.ylabel("Apparent Magnitude")

pl.legend( (p1,p2), ("Data","Interpolation"), loc=2, frameon=False)

#pl.show()

pl.savefig('Lagrange.pdf')

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl

# get and slice data
data = np.loadtxt('cepheids.txt', comments= '#')
t = data[:,0]
m = data[:,1]

# Helper function to compute z
def z(t, t_i, t_i1):
    return (t-t_i)/(t_i1 - t_i)

# Helper function to compute Phi_0
def phi0(z):
    return 2*z*z*z - 3*z*z + 1

# Helper function to compute Phi_1
def phi1(z):
    return z*z*z - 2*z*z + z

# Helper function to compute a list of forward differences
def for_dif(t, m):
    return (m[1:] - m[:-1])/(t[1:] - t[:-1])

# Function to compute piecewise cubic Hermite interpolation
def hermite(x, t):
    y = [m[0]]
    primes = for_dif(t, m)
    for i, point in enumerate(t):
        use = [val for val in x if (t[i] < val <= t[i+1])] #select interval
        if (i+1) == (len(t)-1): #for when we reach final interval
            for num in use:
                z_val = z(num, t[i], t[i+1])
                first = m[i]*phi0(z_val) + m[i+1]*phi0(1-z_val)
                second = primes[i]*(t[i+1]-t[i])*phi1(z_val)
                third = primes[i]*(t[i+1] - t[i])*phi1(1-z_val)
                value = first + second - third
                y.append(value)
        else:
            for num in use:
                z_val = z(num, t[i], t[i+1])
                first = m[i]*phi0(z_val) + m[i+1]*phi0(1-z_val)
                second = primes[i]*(t[i+1]-t[i])*phi1(z_val)
                third = primes[i+1]*(t[i+1] - t[i])*phi1(1-z_val)
                value = first + second - third
                y.append(value)
    return y

x = np.arange(0.0, 1.0, 0.001)
y = hermite(x, t)

# plot data and interpolation
p1, = pl.plot(t, m, 'ro', linewidth=2)
p2, = pl.plot(x, y, "b", linewidth=2)

# set x and y ranges
pl.xlim(min(t),max(t))
pl.ylim(0.0,max(y)*1.05)

# label the axes
pl.xlabel("Time")
pl.ylabel("Apparent Magnitude")

pl.legend( (p1,p2), ("Data","Hermite"), loc=2, frameon=False)

#pl.show()

pl.savefig('Hermite.pdf')

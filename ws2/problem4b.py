#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl

# get and slice data
data = np.loadtxt('cepheids.txt', comments= '#')
t = data[:,0]
m = data[:,1]

# Function to compute a linear piecewise interpolation
def linear(x, t):
    y = [m[0]]
    for i, point in enumerate(t):
        use = [val for val in x if (t[i] < val <= t[i+1])] #list comprehension to select interval
        for num in use:
            p = m[i] + (m[i+1] - m[i])*(num - t[i])/(t[i+1] - t[i]) #compute value
            y.append(p)
    return y

# Function to compute a quadratic piecewise interpolation
def quad(x, t):
    y = [m[0]]
    for i, point in enumerate(t):
        use = [val for val in x if (t[i] < val <= t[i+1])] #list comprehension to select interval
        # need to change point selection depending on the interval
        if (i+1) == (len(t)-1):
            for num in use:
                first = m[i-1]*(num - t[i])*(num - t[i+1])/((t[i-1] - t[i])*(t[i-1]-t[i+1]))
                second = m[i]*(num - t[i-1])*(num - t[i+1])/((t[i] - t[i-1])*(t[i]-t[i+1]))
                third = m[i+1]*(num - t[i-1])*(num - t[i])/((t[i+1] - t[i-1])*(t[i+1]-t[i]))
                y.append(first + second + third)
        else:
            for num in use:
                first = m[i]*(num - t[i+1])*(num - t[i+2])/((t[i] - t[i+1])*(t[i]-t[i+2]))
                second = m[i+1]*(num - t[i])*(num - t[i+2])/((t[i+1] - t[i])*(t[i+1]-t[i+2]))
                third = m[i+2]*(num - t[i])*(num - t[i+1])/((t[i+2] - t[i])*(t[i+2]-t[i+1]))
                y.append(first + second + third)
    return y

x = np.arange(0.0, 1.0, 0.001)
y = linear(x, t)
y2 = quad(x, t)

# plot data and interpolation
p1, = pl.plot(t, m, 'ro', linewidth=2)
p2, = pl.plot(x, y, "b", linewidth=2)
p3, = pl.plot(x, y2, "g", linewidth=2)

# set x and y ranges
pl.xlim(min(t),max(t))
pl.ylim(0.0,max(y2)*1.05)

# label the axes
pl.xlabel("Time")
pl.ylabel("Apparent Magnitude")

pl.legend( (p1,p2,p3), ("Data","Linear", "Quadratic"), loc=2, frameon=False)

#pl.show()

pl.savefig('Piecewise.pdf')

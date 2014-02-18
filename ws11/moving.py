#!/usr/bin/env python

import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp

def analytic(x, x0, v, t, sigma):
    in_exp = -((x-x0-v*t)**2)/(2*sigma*sigma)
    return np.exp(in_exp)

# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = np.linspace(0,100,num=1000)       # CHANGE I'VE   MADE
# parameters
dx = x[1]-x[0]
v = 0.1

n = len(x)
y = np.zeros(n)
cfl = 1.0
dt = cfl*dx/np.abs(v)           # CHANGE I'VE MADE
t = 0.0

# for initial data
sigma = np.sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = analytic(x, x0, v, t, sigma)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x, y, 'r-')
mpl.show()

yold2 = y
yold = y
ntmax = 2000

for it in range(ntmax):
    t = t + dt

    # Save old data
    yold2 = yold
    yold = y

    # Calculate new data
    y = analytic(x, x0, v, t, sigma)

    mpl.clf()
    mpl.plot(x, y, 'r-')
    mpl.draw()

mpl.show()

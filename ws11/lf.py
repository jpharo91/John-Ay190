#!/usr/bin/env python

import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp

def apply_bcs(x,y):
    # apply boundary conditions
    # you need to fill in code
    y[1] = y[0]
    y[-1] = y[-2]
    return y

def analytic(x, x0, v, t, sigma):
    in_exp = -((x-x0-v*t)**2)/(2*sigma*sigma)
    return np.exp(in_exp)

# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = np.linspace(0,100,num=1000)
# parameters
dx = x[1]-x[0]
v = 0.1

n = len(x)
y = np.zeros(n)
cfl = 1.0
dt = 0.1
t = 0.0

# for initial data
sigma = np.sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = analytic(x,x0,v,t,sigma)

times = [t]
errors = [0]

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.subplot(211)
p1, =mpl.plot(x,y,'x-') # numerical data
p2, =mpl.plot(x,analytic(x,x0,v,t,sigma),'r-') # analytic data
mpl.xlabel("x")
mpl.ylabel("y")
mpl.subplot(212)
mpl.plot(times, errors, 'g-')
mpl.xlabel("Time")
mpl.ylabel("Absolute Difference")
mpl.show()

yold2 = y
yold = y
ntmax = 500

for it in range(ntmax):
    t = t + dt
    times.append(t)
    # save previous and previous previous data
    yold2 = yold
    yold = y

    # get new data; ideally just call a function
    #y = ????
    du = []
    du2 = []
    du = np.append(0, yold[2:]-yold[:-2])
    du = np.append(du,0)
    du2 = np.append(0, yold[2:]+yold[:-2])
    du2 = np.append(du2,0)
    y = 0.5*du2 - (v*dt/(2.0*dx))*du

    # after update, apply boundary conditions
    y = apply_bcs(x,y) 

    # get analytic result for time t
    yana = analytic(x, x0, v, t, sigma)
    # compute error estimage
    err = max(np.abs(y-yana))
    errors.append(err)
    print "it = ",it,err
    mpl.clf()
    mpl.subplot(211)
    # plot numerical result
    p1, =mpl.plot(x,y,'x-')
    # plot analytic results
    p2, =mpl.plot(x,yana,'r-')
    mpl.xlabel("x")
    mpl.ylabel("y")
    mpl.subplot(212)
    mpl.plot(times, errors, 'g-')
    mpl.xlabel("Time")
    mpl.ylabel("Absolute Difference")
    mpl.draw()


mpl.savefig("LF4.pdf")

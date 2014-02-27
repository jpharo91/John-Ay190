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


def analytic(x, v, t):
    return (1.0/8.0)*np.sin(2.0*np.pi*(x - v*t)/100.0)

# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = np.linspace(0, 100, num=1000)
# parameters
dx = x[1]-x[0]
u = 0.1

n = len(x)
y = np.zeros(n)
cfl = 1.0
t = 0.0

#set up initial conditions
y = analytic(x, u, t)

times = [t]
errors = [0]
steps = []

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'x-') # numerical data
mpl.plot(x,analytic(x, u, t),'r-') # analytic data
mpl.show()

u = y
ntmax = 200
for it in range(ntmax):
    dt = cfl*dx/max(np.abs(u))
    steps.append(dt)
    t = t + dt
    times.append(t)
    # save previous and previous previous data
    u = y

    # get new data; ideally just call a function
    #y = ????
    y = []
    yana = []
    diff = []
    for i,val in enumerate(u):
        yana_i = analytic(x[i], val, t)
        if val < 0:
            if i < (len(u)-1):
                y_i = val - (val*dt/dx)*(u[i+1]-u[i]) #downwind
            else:
                y_i = 0
        else:
            if i > 0:
                y_i = val - (val*dt/dx)*(u[i]-u[i-1]) #upwind
            else:
                y_i = 0
        diff.append(y_i-yana_i)
        y.append(y_i)
        yana.append(yana_i)

    # after update, apply boundary conditions
    # apply_bcs(x,y) 
    y = apply_bcs(x,y)

    # get analytic result for time t
   # yana = analytic(x, u, t)
    # compute error estimage
    abs_diff = map(np.abs,diff)
    err = max(abs_diff)
    errors.append(err)
    
    print "it = ",it,err
    mpl.clf()
    mpl.subplot(211)
    # plot numerical result
    p1, =mpl.plot(x, y, 'x-')
    # plot analytic results
    mpl.plot(x,yana,'r-')
    mpl.xlabel("x")
    mpl.ylabel("y")
    mpl.subplot(212)
    mpl.plot(times[1:], steps, 'g-')
    mpl.xlabel("Time")
    mpl.ylabel("Time Step")
    mpl.draw()


mpl.savefig("t=200.pdf")



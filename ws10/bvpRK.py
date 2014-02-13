#!/user/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# set up grid
xmin = 0.0
xmax = 1.0
npoints = 1000
# set up grid
x = np.linspace(xmin, xmax, npoints)
# dx based on x[1] and x[0]
dx = x[1] - x[0]

# boundary values
A = 0.0 # inner boundary
B = 0.1 # outer boundary

def calc_rhs(u,xx):
    # rhs routine
    # rhs[0] is rhs for y
    # rhs[1] is rhs for u
    rhs = np.zeros(2)
    rhs[0] = u 
    rhs[1] = 12*xx - 4

    return rhs

def integrate_FE(z,x):
    # forward-Euler integrator
    
    # make an array for all points
    # entry 0 contains y
    # entry 1 contains y'
    yy = np.zeros((npoints,2))

    yy[0,0] = A # boundary value A for y at x=0
    yy[0,1] = z # guessed boundary value for y' at x=0

    for i in range(npoints-1):
        yy[i+1,:] = yy[i,:] + dx*calc_rhs(yy[i,1],x[i])

    return yy

def integrate_RK2(z, x):
    # Runge-Kutta 2 integrator

    # make an array for all points
    # entry 0 contains y
    # entry 1 contains y'
    yy = np.zeros((npoints, 2))

    yy[0,0] = A # boundary value A for y at x=0
    yy[0,1] = z # guessed boundary value for y' at x=0

    for i in range(npoints-1):
        k1 = dx*calc_rhs(yy[i,1],x[i])
        k2 = dx*calc_rhs(yy[i,1]  + 0.5*k1[0], x[i] + 0.5*dx)
        yy[i+1,:] = yy[i,:] + k2

    return yy

# get initial guess for derivative
z0 = -1100000.0
z1 = 10000000.0
yy0 = integrate_RK2(z0,x)
yy1 = integrate_RK2(z1,x)
phi0 = yy0[npoints-1,0] - B
phi1 = yy1[npoints-1,0] - B
dphidz = (phi1 - phi0) / (z1 - z0) # dphi/dz

i = 0
itmax = 100
err = 1.0e99
criterion = 1.0e-12

z0 = z1
phi0 = phi1
while (err > 1.0e-12 and i < itmax):
    z1 = z0 - phi0 / dphidz # secant update
    yy = integrate_RK2(z1,x)
    phi1 = yy[npoints-1,0] - B
    dphidz = (phi1 - phi0) / (z1 - z0) # dphi/dz numerical
    err = np.abs(phi1) # your error measure
    z0 = z1
    phi0 = phi1
    i = i+1

    print i,z1,phi1

# Function for analytic solution
def f(x):
    return 2.0*x**3 - 2*x**2 + 0.1*x

funcvals = []
difs = []

# Calculate absolute difference
for num in x:
    funcvals.append(f(num))

difs = np.abs(funcvals - yy[:,0])

plt.figure(1)
plt.subplot(211)
p1, =plt.plot(x,yy[:,0],"r-")
p2, =plt.plot(x,2.0*x**3 - 2*x**2 + 0.1*x,"k-")
plt.xlabel("x")
plt.ylabel("y")
plt.legend((p1,p2),("Numerical","Analytic"),loc=3,frameon=False)
plt.subplot(212)
plt.plot(x,difs,"g-")
plt.xlabel("x")
plt.ylabel("Absolute Difference")
plt.savefig("BVP-RK2.pdf")

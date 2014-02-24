#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from scipy import interpolate

# global constants
grav = 6.67e-8
msun = 1.99e33

# load data table
data = np.loadtxt('presupernova.dat')

# slice data
index = data[:,0]
mass = data[:,1]
radius = data[:,2]
temp = data[:,3]
density = data[:,4]
vel = data[:,5]
FE = data[:,6]
zeroes = data[:,7]

# set up grid
min_rad =  min(radius)
max_rad = 10**9
npoints = 10000

def index(list, max):
    for i in range(len(list)):
        if list[i] > max:
            return i
    return 0

n_max = index(radius, 10e8)
new_dense = density[:n_max+1]
new_rad = radius[:n_max+1]
#new_mass = mass[:n_max+1]

def rhs(r, phi, z, rho):
    if r == 0.:
        r = 1e-16

    dz = 4.0*np.pi*grav*rho - 2.0*z/r
    dphi = z
    dm = 4.0*np.pi*(r**2)*rho

    return np.array([dz, dphi, dm])

def integrate_FE(r, dr, phi, z, m, rho):
    old = np.zeros(3)
    new = np.zeros(3)

    old[0] = z
    old[1] = phi
    old[2] = m

    new = old + dr*rhs(r, phi, z, rho)
    return new

new_z = np.zeros(npoints)
new_phi = np.zeros(npoints)
new_mass = np.zeros(npoints)

new_z[0] = 0
new_phi[0] = 0
new_mass[0] = mass[0]

rad = np.linspace(new_rad[0], max_rad, num=npoints)
dr = rad[1]-rad[0]
rhos = interpolate.interp1d(new_rad, new_dense, kind='cubic')
densities = rhos(rad)

for i in range(npoints-1):
    new = integrate_FE(rad[i], dr, new_phi[i], new_z[i], new_mass[i], densities[i])
    new_z[i+1] = new[0]
    new_phi[i+1] = new[1]
    new_mass[i+1] = new[2]

print new_phi

const = -grav*new_mass[-1]/rad[-1]
new_phi += const-new_phi[-1]
print new_phi

def phi2(rho, r):
    return (2.0*np.pi/3.0)*grav*rho*(r*r - 3.0*max_rad*max_rad)

analytic = phi2(densities, rad)
print analytic

ax = pl.gca()
ax.set_xscale('Log')
#ax.set_yscale('Log')
p1, = pl.plot(rad, new_phi, 'bo')
p2, = pl.plot(rad, analytic, 'r')
pl.show()

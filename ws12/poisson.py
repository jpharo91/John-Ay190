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

# Function to find index of radius cutoff
def index(list, max):
    for i in range(len(list)):
        if list[i] > max:
            return i
    return 0

# Slice data at radius cutoff
n_max = index(radius, 10e8)
new_dense = density[:n_max+1]
new_rad = radius[:n_max+1]

# Right-hand-side function for calculating the potential
def rhs(r, phi, z, rho):
    if r == 0.:
        r = 1e-16

    dz = 4.0*np.pi*grav*rho - 2.0*z/r
    dphi = z
    dm = 4.0*np.pi*(r**2)*rho

    return np.array([dz, dphi, dm])

# Function for finding the potential via the Forward-Euler method
def integrate_FE(r, dr, phi, z, m, rho):
    old = np.zeros(3)
    new = np.zeros(3)

    old[0] = z
    old[1] = phi
    old[2] = m

    new = old + dr*rhs(r, phi, z, rho)
    return new

# Set up arrays
new_z = np.zeros(npoints)
new_phi = np.zeros(npoints)
new_mass = np.zeros(npoints)

# Initial conditions
new_z[0] = 0
new_phi[0] = 0
new_mass[0] = mass[0]

# Generate radius and density grids via interpolation
rad = np.linspace(new_rad[0], max_rad, num=npoints)
dr = rad[1]-rad[0]
rhos = interpolate.interp1d(new_rad, new_dense, kind='cubic')
rho = 10**5 # Select a constant density for use in the homogeneous sphere
densities = rho * np.ones(npoints)

# Use FE to calculate phi for homogeneous sphere
for i in range(npoints-1):
    new = integrate_FE(rad[i], dr, new_phi[i], new_z[i], new_mass[i], densities[i])
    new_z[i+1] = new[0]
    new_phi[i+1] = new[1]
    new_mass[i+1] = new[2]

# Adjust constant
const = -grav*new_mass[-1]/rad[-1]
new_phi += const-new_phi[-1]

# Function for Equation 3
def ana_phi(rho, r):
    return (2.0*np.pi/3.0)*grav*rho*(r*r - 3.0*max_rad*max_rad)

analytic = ana_phi(densities, rad)

rel_err = np.abs((analytic - new_phi)/analytic)

# Repeat the above steps, but with twice as many grid points
rad2 = np.linspace(new_rad[0], max_rad, num=2*npoints)
dr2 = rad2[1]-rad2[0]
densities2 = rho * np.ones(2*npoints)

z2 = np.zeros(2*npoints)
phi2 = np.zeros(2*npoints)
mass2 = np.zeros(2*npoints)

z2[0] = 0
phi2[0] = 0
mass2[0] = mass[0]

for i in range(2*npoints - 1):
    new3 = integrate_FE(rad2[i], dr2, phi2[i], z2[i], mass2[i], densities2[i])
    z2[i+1] = new3[0]
    phi2[i+1] = new3[1]
    mass2[i+1] = new3[2]

constant = -grav*mass2[-1]/rad2[-1]
phi2 += constant - phi2[-1]

analytic2 = ana_phi(densities2, rad2)

rel_err2 = np.abs((analytic2 - phi2)/analytic2)

rate = rel_err2[0]/rel_err[0]
print "Convergence Rate: " + str(rate)

# Apply code to the presupernova model
model_rho = rhos(rad)
model_z = np.zeros(npoints)
model_phi = np.zeros(npoints)
model_mass = np.zeros(npoints)

model_z[0] = 0
model_phi[0] = 0
model_mass[0] = mass[0]

for i in range(npoints-1):
    new2 = integrate_FE(rad[i], dr, model_phi[i], model_z[i], model_mass[i], model_rho[i])
    model_z[i+1] = new2[0]
    model_phi[i+1] = new2[1]
    model_mass[i+1] = new2[2]

const2 = -grav*model_mass[-1]/rad[-1]
model_phi += const2 - model_phi[-1]

ax = pl.gca()
ax.set_xscale('Log')
#ax.set_yscale('Log')
#pl.figure(1)
#pl.subplot(211)
#p1, = pl.plot(rad, new_phi, 'b--', lw=2)
#p2, = pl.plot(rad, analytic, 'r', lw=2)
#pl.ylabel("Potential (erg/g)")
#pl.legend((p1,p2), ("Numerical","Analytic"), loc=2)
#pl.subplot(212)
#p3, = pl.plot(rad, rel_err, 'g', lw=2)
#pl.xlabel("Radius (cm)")
#pl.ylabel("Relative Error")
#pl.subplot(222)
p4, = pl.plot(rad, model_phi, 'k--', lw=2)
pl.xlabel("Radius (cm)")
pl.ylabel("Potential (erg/g)")
pl.savefig("Presupernova.pdf")

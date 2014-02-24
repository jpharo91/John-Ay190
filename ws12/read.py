#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl

data = np.genfromtxt('presupernova.dat')
num = data[:,0]
mass = data[:,1] #increase, goes to 10e34, so must be enclosed mass
third = data[:,2] #increase, goes to ~500 solar radii, so must be radius
fourth = data[:,3] #must be temp b/c only decreasing one apart from density
fifth = data[:,4] #must be density b/c only monotonically decreasing
sixth = data[:,5] # radial velocity profile
seventh = data[:,6] # number of electrons per baryon, only important in core
eighth = data[:,7] #just zeroes?

ax = pl.gca()
#p1, = pl.plot(num, second, 'r-')
#p2, = pl.plot(num, third, 'b-')
#p3, = pl.plot(third, fourth, 'g-')
#p4, = pl.plot(num, fifth, 'k-')
#p5, = pl.plot(third, sixth)
#p6, = pl.plot(third, seventh)
#p7, = pl.plot(third, seventh)
#ax.set_xscale('Log')
#ax.set_yscale('Log')
#pl.xlabel("Radius (cm)")
#pl.ylabel("???")
#pl.savefig("EighthColumn.pdf")

max_ind = 0
max_rad = 10e8

for i,rad in enumerate(third):
    if rad >= max_rad:
        max_ind = i
        break

core_rad = third[:max_ind]
core_rho = fifth[:max_ind]

ax.set_xscale('Log')
ax.set_yscale('Log')
pl.plot(core_rad, core_rho)
pl.xlabel("Radius (cm)")
pl.ylabel("Density (g/cm^3")
pl.savefig("Density-Radius.pdf")

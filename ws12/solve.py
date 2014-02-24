#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from scipy import interpolate

data = np.loadtxt('presupernova.dat')

index = data[:,0]
mass = data[:,1]
radius = data[:,2]
temp = data[:,3]
density = data[:,4]
vel = data[:,5]
FE = data[:,6]
zeroes = data[:,7]

min_rad =  min(radius)
max_rad = 10e8
npoints = 10000.0
dr = (max_rad - min_rad)/npoints

def index(list, max):
    for i in range(len(list)):
        if list[i] > max:
            return i
    return 0

n_max = index(radius, 10e8)
new_dense = density[:n_max+1]
new_rad = radius[:n_max+1]

x = np.arange(min_rad, max_rad, dr)
y = interpolate.interp1d(new_rad, new_dense, kind='cubic')

ax = pl.gca()
ax.set_xscale('Log')
ax.set_yscale('Log')
#p1, = pl.plot(new_rad, new_dense, 'kx')
pl.plot(x, y(x), 'b')
pl.xlim(min(x)*0.95,max(x)*1.05)
pl.xlabel("Radius (cm)")
pl.ylabel("Interpolated Density (g/cm^3)")
pl.savefig("Interpolation.pdf")

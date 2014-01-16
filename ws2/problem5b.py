#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import interp1d

# get and slice data
data = np.loadtxt('cepheids.txt', comments= '#')
t = data[:,0]
m = data[:,1]

x = np.arange(0.0, 1.0, 0.001)
y = interp1d(t, m, kind='cubic')

# plot data and interpolation
p1, = pl.plot(t, m, 'ro', linewidth=2)
p2, = pl.plot(x, y(x), "b", linewidth=2)

# set x and y ranges
pl.xlim(min(t),max(t))
#pl.ylim(0.0,max(y)*1.05)

# label the axes
pl.xlabel("Time")
pl.ylabel("Apparent Magnitude")

pl.legend( (p1,p2), ("Data","Cubic Spline"), loc=2, frameon=False)

#pl.show()

pl.savefig('Spline.pdf')

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl

# read in data
data = np.genfromtxt("modified_m_sigma_table.dat")

# separate data
redshift = data[:,0]
velocity = data[:,1]
sigma_uncertainty = data[:,2]
h_alpha = data[:,3]
h_error = data[:,4]
luminosity = data[:,5]
l_error = data[:,6]
BHM = data[:,7]
BHM_upper = data[:,8]
BHM_lower = data[:,9] # formal uncertainty

# Take log of velocity dispersion
log_velocity = np.log10(velocity)

# Make basic plot of data
p1, = pl.plot(log_velocity, BHM, "bo", linewidth=2)

pl.xlabel("Log Velocity Dispersion (km/s)")
pl.ylabel("Log Black Hole Mass (in solar masses)")

#pl.show()

pl.savefig("PartA.pdf")

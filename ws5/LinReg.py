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

# Take log of velocity dispersion and error
log_velocity = np.log10(velocity)
log_sig_vel = np.log10(sigma_uncertainty)

# list of constant errors
const_sig = np.ones(len(BHM))

# Helper function to compute the values of the fit, taking only one list of errors
def linreg2(x, y, sig):
    S = sum(1/(sig*sig))
    sumx = sum(x/(sig*sig))
    sumx2 = sum((x*x)/(sig*sig))
    sumy = sum(y/(sig*sig))
    sumy2 = sum((y*y)/(sig*sig))
    sumxy = sum((x*y)/(sig*sig))
    a1 = (sumy*sumx2 - sumx*sumxy)/(S*sumx2 - sumx*sumx)
    sig_a1 = np.sqrt(sumx2/(S*sumx2 - (sumx*sumx)))
    a2 = (S*sumxy - sumy*sumx)/(S*sumx2 - sumx*sumx)
    sig_a2 = np.sqrt(S/(S*sumx2 - (sumx*sumx)))
    chi2 = sum((a1 + a2*x - y)**2/(sig*sig))
    return (a1, sig_a1, a2, sig_a2, chi2)

# Function to calculate linear regression including x and y errors
def linreg(x, y, sigx, sigy):
    (tempA1, temp_sig1, tempA2, temp_sig2, tempChi2) = linreg2(x, y, sigy)
    sig_ex = tempA2*sigx
    sig = sigy*sigy + sig_ex*sig_ex
    return linreg2(x, y, sig)

# Using constant sigma, so errors cancel. Chi^2 values is meaningless
#(a1, sig_a1, a2, sig_a2, chi2) = linreg2(log_velocity, BHM, const_sig)
#print (a1, sig_a1, a2, sig_a2, chi2)

# Variable sigma results
(a1, sig_a1, a2, sig_a2, chi2) = linreg(log_velocity, BHM, log_sig_vel, BHM_lower)
print (a1, sig_a1, a2, sig_a2, chi2)

# Apply fit in order to plot regression
xs = np.arange(min(log_velocity)-5, max(log_velocity)+5, .01)
ys = a1 + a2*xs

p1, = pl.plot(log_velocity, BHM, "bo", linewidth=2)
p2, = pl.plot(xs, ys, "r", linewidth=2)

pl.xlabel("Log Velocity Dispersion (km/s)")
pl.ylabel("Log Black Hole Mass (in solar masses)")

# Add error bars
ax = pl.gca()
ax.errorbar(log_velocity, BHM, log_sig_vel, BHM_lower, color="b", linestyle="None", marker="None")

# set x and y ranges
#pl.xlim(min(log_velocity)*.9,max(log_velocity)*1.1)
#pl.ylim(min(BHM*.9),max(BHM*1.1))
pl.xlim(1.4,2.6)
pl.ylim(4.5,9.5)

pl.legend( (p1,p2), ("Data","Linear Fit"), loc=2, frameon=False )

#pl.show()

pl.savefig("PartC.pdf")

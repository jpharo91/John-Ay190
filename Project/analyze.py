#!/usr/bin/env python
    
"""
John Pharo and Cutter Coryell
Caltech Ay 190 Winter 2014 with Christian Ott
Final Project

Analyzes data produced by simulate.py.
"""

import numpy as np
import matplotlib.pyplot as pl
import cPickle as pickle
import subprocess
import plot_defaults


####################################
# General Constants (CGS units)
ggrav = 6.67e-8
c = 3.0e10
msun = 1.99e33
rsun = 6.96e10
direcs = [0, 1, 2] # indices corresponding to the three spatial directions

####################################
# Parameters

output_dir = "example" # name of directory containing simulation output data
start_iteration = 0 # the first iteration to load
n_points = 1e2 # number of iterations in the simulation
input_every = 1 # inputs one data file for every `input_every` data files
m1 = 1.4*msun # mass of body 1
m2 = 1.4*msun # mass of body 2

####################################
# Load Data

n_loaded = (n_points - start_iteration) / input_every # number of data files to load

# array of time at each time step
times = np.zeros(n_loaded)

# array of true anomaly at each time step
phis = np.zeros(n_loaded)

# array of semi-major axis at each time step
axs = np.zeros(n_loaded)

# array of eccentricity at each time step
eccs = np.zeros(n_loaded)

# array of separation vector at each time step
xs = np.zeros((3, n_loaded))

# array of reduced quadrupole tensor at each time step
I_bars = np.zeros((3, 3, n_loaded))

# array of plus-polarized gravitational wave strain at each time step
h_pluses = np.zeros(n_loaded)

# array of cross-polarized gravitational wave strain at each time step
h_crosses = np.zeros(n_loaded)

for it in range(int(n_loaded)):
    f = open("{}/iteration-{}.pickle".format(output_dir, start_iteration + it*input_every), 'rb')
    data = pickle.load(f)
    f.close()
    
    times[it], phis[it], axs[it] = data['t'], data['phi'], data['a']
    eccs[it], xs[:,it] = data['e'], data['x']
    h_pluses[it], h_crosses[it] = data['h_plus'], data['h_cross']
    try:
        I_bars[:,:,it] = data['I_bar2']
    except KeyError: # backwards compatibility with old format of output data
        I_bars[:,:,it] = data['I_bar']
    try:
        m1, m2 = data['m1'], data['m2']
    except KeyError:
        pass # use the parameter values defined in this file

####################################
# Analyze Data

M = m1 + m2 # total mass
mu = (m1 * m2) / M # reduced mass

# Calculate z-component of orbital angular momentum
Ls = np.sqrt(ggrav * M * axs * (1 - eccs**2))

# Calculate total mechanical energy of the system
Es = 0.5 * (ggrav * M / Ls)**2 * (eccs**2 - 1)

####################################
# Plot Data

fig_dir = output_dir + "_figs"
subprocess.call("mkdir -p " + fig_dir, shell=True)

fig = pl.figure(figsize=(10,8))
fig.subplots_adjust(left=0.2)
fig.subplots_adjust(bottom=0.16)
fig.subplots_adjust(top=0.85)
fig.subplots_adjust(right=0.85)

pl.clf()
pl.plot(times, np.linalg.norm(xs, axis=0), lw=8)
pl.xlabel(r"Time \(t\) [s]")
pl.ylabel(r"Radius \(r\) [cm]")
pl.xlim(times[0], times[-1])
pl.savefig(fig_dir + "/radius.pdf")

pl.clf()
pl.plot(times, phis, lw=8)
pl.xlabel(r"Time \(t\) [s]")
pl.ylabel(r"True anomaly \(\phi\) [radians]")
pl.xlim(times[0], times[-1])
pl.savefig(fig_dir + "/true-anomaly.pdf")

pl.clf()
h_plus, = pl.plot(times[1:], h_pluses[1:] * ggrav / c**4, lw=8)
h_cross, = pl.plot(times[1:], h_crosses[1:] * ggrav / c**4, 'r', lw=8)
pl.legend( (h_plus, h_cross), (r"$h_+$", r"$h_{\times}$"), frameon=False, loc="lower left")
pl.xlabel(r"Time \(t\) [s]")
pl.ylabel(r"Gravitational Wave Strain")
pl.xlim(times[0], times[-1])
pl.savefig(fig_dir + "/gw-strain.pdf")

pl.clf()
pl.plot(times, axs, lw=8)
pl.xlabel(r"Time \(t\) [s]")
pl.ylabel(r"Semi-major axis \(a\) [cm]")
pl.xlim(times[0], times[-1])
pl.savefig(fig_dir + "/semi-major-axis.pdf")

pl.clf()
pl.plot(times, eccs, lw=8)
pl.xlabel(r"Time \(t\) [s]")
pl.ylabel(r"Eccentricity \(e\)")
pl.xlim(times[0], times[-1])
pl.savefig(fig_dir + "/eccentricity.pdf")

pl.clf()
pl.plot(times, Es, lw=8)
pl.xlabel(r"Time \(t\) [s]")
pl.ylabel(r"Energy \(E\) [erg]")
pl.xlim(times[0], times[-1])
pl.savefig(fig_dir + "/energy.pdf")

pl.clf()
pl.plot(times, Ls, lw=8)
pl.xlabel(r"Time \(t\) [s]")
pl.ylabel(r"Angular momentum \(L_z\) [erg \(\cdot\) s]")
pl.xlim(times[0], times[-1])
pl.savefig(fig_dir + "/angular-momentum.pdf")

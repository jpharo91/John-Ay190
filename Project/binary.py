#!/usr/bin/env python

import numpy as np

# general constants
ggrav = 6.67e-8
c = 3.0e10
msun = 1.99e33


####################################
# helper functions

# Function to calculate the first time derivative of `quantity`
def first_time_der(quantity):

# Function to calculate the secoond time derivative of `quantity`
def second_time_der(quantity):

# Function to calculate the third time derivative of `quantity`
def third_time_der(quantity):

# Function to calculate new orbital separation vector
def orbit():

# Function to calculate current reduced quadrupole tensor
def quadrupole():

# Function to calculate the righthand-side of the energy evolution equation
def energy_RHS():

# Function to calculate the righthand-side of the angular momentum evolution equation
def momentum_RHS():

# Function to calculate the righthand-side of the energy evolution equation
def integrate_energy():

# Function to calculate the righthand-side of the angular momentum evolution equation
def integrate_momentum():

# Function to calculate gravitational wave
def wave():


####################################
# main code body

# parameters
npoints = 10000
t = 0.0
t_final = 100.0
m1 = 1.4*msun
m2 = 1.0*msun

# set up initial conditions
M = m1 + m2
mu = (m1 * m2) / M

times = np.linspace(t, t_final, num=npoints)
dt = times[1] - times[0]
x = [1.0, 1.0, 1.0]
x_mag = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
x_min = 2e6

# total system energy
energy = np.zeros(npoints)
# total system angular momentum
momentum = np.zeros(npoints)
# system reduced quadrupole tensor (sets up a 3x3 dictionary indexed like quad[i][j])
quad = {1:{}, 2:{}, 3:{}}
for i in range(1,4):
    quad[i] = {}
    for j in range(1,4):
        quad[i][j] = np.zeros(npoints)


# main loop
for it, t in enumerate(times):
    energy[it] = integrate_energy()
    momentum[it] = integrate_momentum()
    x = orbit()
    x_mag = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
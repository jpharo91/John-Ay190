#!/usr/bin/env python

import numpy as np

# general constants
ggrav = 6.67e-8
c = 3.0e9
msun = 1.99e33


####################################
# helper functions

# Function to calculate array of reduced quadrupole tensor
def quadrupole():

# Function to calculate time derivatives of the reduced quadrupole tensor
def dquadrupole():

# Function to calculate the energy loss rate
def energy_loss():

# Function to caculate the momentum loss rate
def momentum_loss():

# Function to calculate timestep
def timestep():

# Function to calculate new energy
def energy():

# Function to calculate new energy
def momentum():

# Function to calculate new orbital separation vector
def orbit():

# Function to calculate gravitational wave
def wave():


####################################
# main code body

# set up initial conditions
m1 = 1.4*msun
m2 = 1.0*msun
M = m1 + m2
mu = (m1 * m2) / M
t = 0.0
t_final = 5.0

# while loop until end conditions met
# calculate time step, update energy and momentum, update orbits, update time
while (t < t_final):
    dt = timestep
    old_energy = energy
    old_mom = momentum
    energy = old_energy + dt * energy_loss()
    momentum = old_mom + dt * momentum_loss()
    x = orbit()
    t += dt

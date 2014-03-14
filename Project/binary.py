#!/usr/bin/env python

import numpy as np

# general constants
ggrav = 6.67e-8
c = 3.0e10
msun = 1.99e33
frac = 1.0/6.0


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
def energy_RHS(tensor):
    temp = third_time_der(tensor)
    coeff = ggrav / (5.0 * c**5)
    return coeff * temp * temp

def perm_parity(lst):
    '''\
    Given a permutation of the digits 0..N in order as a list, 
    returns its parity (or sign): +1 for even parity; -1 for odd.
    '''
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i,len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity

# Function to calculate the Levi-Civita tensor
def lc(i, j, k):
    if i == j or i == k or j == k:
        return 0
    else:
        return perm_parity([i,j,k])

# Function to calculate the righthand-side of the angular momentum evolution equation
def momentum_RHS(tensor, axis):
    coeff = (2.0 * ggrav)/(5.0 * c**5)
    epsilon = lc(axis, j, k) #Not sure how to make this work
    m_vals = [0,1,2]
    derivs = 0
    for m in m_vals:
        derivs += third_time_der(tensor[j,m]) * third_time_der(tensor[k,m])
    return coeff * epsilon * derivs

# Function to calculate the righthand-side of the energy evolution equation
def integrate_energy(tensor, t, dt, E):
    old = E
    k1 = dt * energyRHS(tensor)
    k2 = dt * energyRHS(t + dt/2.0, tensor + k1/2.0) #Not sure if t is actually necessary
    k3 = dt * energyRHS(t + dt/2.0, tensor + k2/2.0) #Maybe need function to give tensor
    k4 = dt * energyRHS(t + dt, tensor + k3)         # as a function of time

    new = old + frac * (k1 + 2.0*k2 + 2.0*k3 + k4)
    return new

# Function to calculate the righthand-side of the angular momentum evolution equation
def integrate_momentum():
    old = p
    k1 = dt * momentum_RHS()
    k2 = dt * momentum_RHS()
    k3 = dt * momentum_RHS()
    k4 = dt * momentum_RHS()

    new = old + frac * (k1 + 2.0*k2 + 2.0*k3 + k4)
    return new

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

    if x_mag < x_min:
        break
    else:
        energy[it] = integrate_energy()
        momentum[it] = integrate_momentum()
        x = orbit()
        x_mag = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])

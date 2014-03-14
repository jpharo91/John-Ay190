#!/usr/bin/env python

import numpy as np

# general constants
ggrav = 6.67e-8
c = 3.0e10
msun = 1.99e33
frac = 1.0/6.0


####################################
# helper functions


def first_time_der(quantity):
    """
    Function to calculate the first time derivative of `quantity`.
    """


def second_time_der(quantity):
    """
    Function to calculate the second time derivative of `quantity`.
    """

def third_time_der(quantity):
    """
    Function to calculate the third time derivative of `quantity`.
    """


def orbit(E, L, t):
    """
    Function to calculate new orbital separation vector.
    """
    eccentricity = np.sqrt(1.0 + (2.0*E*L*L)/(ggrav*M)**2)
    p = L*L/(ggrav*M)
    r = p/(1 - e*np.cos(phi(t))) #Need to define phi function

def quadrupole():
    """
    Function to calculate current reduced quadrupole tensor.
    """

def energy_RHS(tensor):
    """
    Function to calculate the righthand-side of the energy evolution equation.
    """
    temp = third_time_der(tensor)
    coeff = ggrav / (5. * c**5)
    return coeff * temp * temp

def perm_parity(lst):
    """
    Given a permutation of the digits 0..N in order as a list, 
    returns its parity (or sign): +1 for even parity; -1 for odd.
    """
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i,len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity

def lc(i, j, k):
    """
    Function to calculate the Levi-Civita tensor \( \epsilon_{i,j,k} \).
    """
    if i == j or i == k or j == k:
        return 0
    else:
        return perm_parity([i,j,k])

def momentum_RHS(tensor, axis):
    """
    Function to calculate the righthand-side of the angular momentum evolution equation.
    """
    coeff = (2. * ggrav) / (5. * c**5)
    dirs = [1, 2, 3] # three spatial directions
    total = 0.
    for j in dirs:
        for k in dirs:
            for m in dirs:
                epsilon = lc(axis, j, k)
                total += (coeff * epsilon * second_time_der(tensor[j][m])
                           * third_time_der(tensor[k][m]))
    return coeff * total

def integrate_energy(tensor, t, E):
    """
    Function to calculate the righthand-side of the energy evolution equation
    """
    k1 = energyRHS(quadrupole(t))
    k2 = energyRHS(tensor + k1/2.0) #Not sure if t is actually necessary
    k3 = energyRHS(t + dt/2.0, tensor + k2/2.0) #Maybe need function to give tensor
    k4 = energyRHS(t + dt, tensor + k3)         # as a function of time

    return E + dt * (k1 + 2.*k2 + 2.*k3 + k4) / 6.

def integrate_momentum():
    """
    Function to calculate the righthand-side of the angular momentum evolution equation
    """
    k1 = momentum_RHS()
    k2 = momentum_RHS()
    k3 = momentum_RHS()
    k4 = momentum_RHS()

    return L + dt * (k1 + 2.*k2 + 2.*k3 + k4) / 6.

def wave():
"""
Function to calculate gravitational wave
"""


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

#!/usr/bin/env python

import numpy as np
from scipy import special as sp

# constants
kBT = 20.0 #in MeV
beta = 1.0 / kBT
hbar = 6.582119*10**(-22) #in MeV*s
c = 3.0*10**8 #in m/s

coefficient = 8.0 * np.pi * (kBT)**3 / ((2.0*np.pi*hbar*c)**3)

# Function for the integrand
def f(x):
    return ((x*x)/(1.0 + np.exp(-x)))

# I first wrote this using a for loop, but I abandoned it in the hopes
# that other methods would be faster
#def gl(n):
#    [lg_roots, lg_weights] = sp.l_roots(n,0)
#    sum = 0
#    for i in range(n):
#        sum += lg_weights[i] * f(lg_roots[i])
#        return sum

# Function to calculate Gaussian-Laguerre Quadrature
def gla(n):
    [lag_roots, lag_weights] = sp.l_roots(n,0)
    func = map(f, lag_roots)
    prod = np.dot(lag_weights, func)
    return prod

# Function to calculate number density of electrons
def density(n):
    return coefficient * gla(n)

trial1 = density(10)
trial2 = density(100)
#trial3 = density(1000)  takes a really long time
#trial4 = density(2000)

print trial1
print trial2
#print trial3
#print trial4

coefficient2 = (8.0 * np.pi) / ((2.0 * np.pi * hbar * c)**3)

# Helper function for integrand
def g(E):
    return (E*E)/(np.exp(beta*E) + 1)

# Function to calculate Gaussian-Legendre Quadrature
def gle(n, a, b):
    [leg_roots, leg_weights] = sp.p_roots(n,0)
    # Once again, I abandon the for loop
   # total = 0
   # for i in range(n):
   #     total += leg_weights[i] * g(((b-a)/2)*leg_roots[i] + ((a+b)/2))
   #     return total
    func_vals = g(((b-a)/2)*leg_roots + ((a+b)/2))
    return np.dot(leg_weights, func_vals)

# Function to calculate the electron density by summing the intervals
# created for quadrature
def density2(n):
    total = 0
    for i in range(0, 150, 5):
        total += (((i+5) - i)/2)*gle(n, i, i+5)
    return (coefficient2 * total)*5.0

trial5 = density2(100)
print trial5
trial6 = density2(500)
print trial6

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import mpl_toolkits.plot3d as mpl3d

# general constants (CGS)
ggrav = 6.67e-8
c = 3.0e10
msun = 1.99e33
direcs = [1, 2, 3] # indices corresponding to the three spatial directions


####################################
# helper functions

def kron_delta(i, j):
    """
    Computes the Kronecker delta of `i` and `j`; that is,
    returns 1 if `i` == `j` and 0 otherwise.
    """
    if i == j:
        return 1
    else:
        return 0

def first_time_der(quantity, it):
    """
    Function to calculate the first time derivative of `quantity` at
    iteration `it`.
    """
    return (quantity[it] - quantity[it-1]) / dt

def second_time_der(quantity, it):
    """
    Function to calculate the second time derivative of `quantity` at
    iteration `it`.
    """
    return (quantity[it] - 2*quantity[it-1] + quantity[it-2]) / dt**2

def third_time_der(quantity, it):
    """
    Function to calculate the third time derivative of `quantity` at
    iteration `it`.
    """
    return (quantity[it] - 3*quantity[it-1] + 3*quantity[it-2]
                                              - quantity[it-3]) / dt**3

def eval_e(E, L):
    return np.sqrt(1. + (2. * E * L**2) / (ggrav * M)**2)

def eval_p(L):
    return L**2 / (ggrav * M)

def eval_x(E, L, phi):
    """
    Function to calculate new orbital separation vector.
    """
    e, p = eval_e(E, L), eval_p(L)
    r = p / (1 - e * np.cos(phi))
    return r * np.array([np.cos(phi), np.sin(phi), 0.])


def eval_quad(x):
    """
    Function to calculate new reduced quadrupole tensor.
    """
    I_bar = {1:{}, 2:{}, 3:{}}
    for i in direcs:
        for j in direcs:
            I_bar[i][j] = mu * (x[i-1] * x[j-1] - kron_delta(i,j) * np.sum(x**2) / 3.)
    return I_bar

def phi_RHS(x, L):
    return L / np.sum(x**2)

def integrate_RK4():
    old = new = np.zeros(3)
    old[0], old[1], old[2] = E, L, phi

    k1 = RHS(x, L)

def RHS(x, L):
    return np.array([energy_RHS(), momentum_RHS(), phi_RHS(x, L)])

def integrate_phi(phi, r, L):
    k1 = phi_RHS(r, L)
    k2 = phi_RHS(tensor + k1/2.0) #Not sure if t is actually necessary
    k3 = energyRHS(t + dt/2.0, tensor + k2/2.0) #Maybe need function to give tensor
    k4 = energyRHS(t + dt, tensor + k3)         # as a function of time


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
    direcs = [1, 2, 3] # three spatial directions
    total = 0.
    for j in direcs:
        for k in direcs:
            for m in direcs:
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

# parameters (CGS)
npoints = 10000
t = 0.0 # initial time
t_final = 100.0 # final time
m1 = 1.4*msun # mass of body 1
m2 = 1.0*msun # mass of body 2
init_separation = 1e10 # initial separation
min_separation = 2e6 # simulation terminates if separation falls beneath this

# set up initial conditions
M = m1 + m2
mu = (m1 * m2) / M

times = np.linspace(t, t_final, num=npoints)
dt = times[1] - times[0]
x = init_separation * np.array([1., 0., 0.])

# array of orbital angle at each time step
phis = np.zeros(npoints)
# array total system energy at each time step
energies = np.zeros(npoints)
# array of total system angular momentum at each time step
momenta = np.zeros(npoints)
# 3x3 tensor of arrays of reduced quadrupole components at each time step
quadrupoles = {}
for i in direcs:
    quadrupoles[i] = {}
    for j in direcs:
        quadrupoles[i][j] = np.zeros(npoints)

### main loop ###
for it, t in enumerate(times):

    if np.sum(x**2) < x_min:
        break
    else:
        energy[it] = integrate_energy()
        momentum[it] = integrate_momentum()
        x = orbit()

        # Plot stuff
        plt.clf()
        fig = plt.gcf()
        ax = mpl3d.Axes3D(fig)
        ax.scatter(x[0], x[1], x[2]) #Proably need to make x an np array
        ax.set_xlim((-rmax,rmax)) #Need to define rmax, probably 1.1 times initial x, or something
        ax.set_ylim((-rmax,rmax))
        ax.set_zlim((-rmax,rmax))
        pl.draw()

pl.show()

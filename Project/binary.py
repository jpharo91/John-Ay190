#!/usr/bin/env python

"""
John Pharo and Cutter Coryell
Caltech Ay 190 Winter 2014 with Christian Ott
Final Project

Simulates the orbit of a binary system as it decays due to gravitational wave
emission. 
"""

import numpy as np
import matplotlib.pyplot as pl
import mpl_toolkits.plot3d as mpl3d

####################################
# Constants and Parameters

# general constants (CGS)
ggrav = 6.67e-8
c = 3.0e10
msun = 1.99e33
direcs = [0, 1, 2] # indices corresponding to the three spatial directions

# simulation parameters
npoints = 10000
integrate = integrate_RK4

# system parameters (CGS)
t = 0.0 # initial time
t_final = 100.0 # final time
init_major_axis = 1e10
init_eccentricity = 0.
m1 = 1.4*msun # mass of body 1
m2 = 1.0*msun # mass of body 2
min_separation = 2e6 # simulation terminates if separation falls beneath this


####################################
# Helper Functions

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
    iteration `it`. If not enough data to compute, returns 0.
    """
    try:
        return (quantity[it] - quantity[it-1]) / dt
    except IndexError:
        return 0

def second_time_der(quantity, it):
    """
    Function to calculate the second time derivative of `quantity` at
    iteration `it`.  If not enough data to compute, returns 0.
    """
    try:
        return (quantity[it] - 2*quantity[it-1] + quantity[it-2]) / dt**2
    except IndexError:
        return 0

def third_time_der(quantity, it):
    """
    Function to calculate the third time derivative of `quantity` at
    iteration `it`.  If not enough data to compute, returns 0.
    """
    try:
        return (quantity[it] - 3*quantity[it-1] + 3*quantity[it-2]
                                                  - quantity[it-3]) / dt**3
    except IndexError:
        return 0

def eval_x(phi, a, e):
    """
    Function to calculate new orbital separation vector.
    """
    return eval_r(phi, a, e) * np.array([np.cos(phi), np.sin(phi), 0.])

def eval_quad(x):
    """
    Function to calculate new reduced quadrupole tensor.
    """
    I_bar = np.zeros((3, 3))
    xnormsq = np.sum(x**2)
    for i in direcs:
        for j in direcs:
            I_bar[i,j] = mu * (x[i] * x[j] - kron_delta(i,j) * xnormsq / 3.)
    return I_bar


def integrate_RK4(phi, a, e):
    """
    Integrates the orbital angle phi, the major axis a, and the eccentricity e
    with the 4th-order Runge-Kutta method.
    """
    old = np.array([phi, a, e])
    k1 = RHS(old)
    k2 = RHS(old + 0.5 * k1)
    k3 = RHS(old + 0.5 * k2)
    k4 = RHS(old + k3)
    new = old + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.
    return (new[0], new[1], new[2])

def RHS(quants):
    """
    Provides a length-3 array of the righthand-sides (time-derivatives)
    of phi, a, and e for used in integration schemes. `quants` is a length-3
    sequence with first element phi, second element a, and third element e.
    """
    phi, a, e = quants[0], quants[1], quants[2]
    r = eval_r(phi, a, e)
    return np.array([phi_RHS(r, a, e), a_RHS(a, e), e_RHS(a, e)])

def eval_r(phi, a, e):
    return a * (1 - e**2) / (1 + e * np.cos(phi))

def phi_RHS(r, a, e):
    return np.sqrt(ggrav * M * a * (1 - e**2)) / r**2

def a_RHS(a, e):
    return ((-64./5.) * (ggrav**3 * m1 * m2 * M)
            * (1 + (73./24.) * e**2 + (37./96.) * e**4) 
            / (c**5 * a**3 * (1 - e**2)**3.5))

def e_RHS(a, e):
    return (-304./15.) * e * (ggrav**3 * m1 * m2 * M)
            * (1 + (121./304.) * e**2) / (c**5 * a**4 * (1 - e**2)**2.5)


def wave():
    """
    Function to calculate gravitational wave
    """


####################################
# main code body


# set up initial conditions

M = m1 + m2
mu = (m1 * m2) / M

times = np.linspace(t, t_final, num=npoints)
dt = times[1] - times[0]

# array of orbital angle at each time step
phis = np.zeros(npoints)

# array of major axis at each time step
axs = np.zeros(npoints)
axs[0] = init_major_axis

# array of eccentricity at each time step
eccs = np.zeros(npoints)
eccs[0] = init_eccentricity

# array of separation vector at each time step
xs = np.zeros((3, npoints))
xs[:,0] = eval_x(phis[0], axs[0], eccs[0])

# array of reduced quadrupole tensor at each time step
quads = np.zeros((3, 3, npoints))

# array of plus-polarized gravitational wave strain at each time step
hpluses = np.zeros(npoints)

# array of cross-polarized gravitational wave strain at each time step
hcrosses = np.zeros(npoints)


### main loop ###
for it, t in enumerate(times[:,-1]):

    # Update quantities
    phis[it+1], axs[it+1], eccs[it+1] = integrate(phis[it], axs[it], eccs[it])
    xs[:,it+1] = eval_x(phis[it+1], axs[it+1], eccs[it+1])
    quads[:,:,it+1] = eval_quad(xs[it+1])
    r = np.linalg.norm(xs[:,it+1])
    hpluses[it+1] = second_time_der(quads[0,0,:] - quads[1,1,:], it+1) / r
    hcrosses[it+1] = 2 * second_time_der(quads[0,1,:]) / r

    # Plot stuff
    plt.clf()
    fig = plt.gcf()
    ax = mpl3d.Axes3D(fig)
    ax.scatter(x[0], x[1], x[2]) #Proably need to make x an np array
    ax.set_xlim((-rmax,rmax)) #Need to define rmax, probably 1.1 times initial x, or something
    ax.set_ylim((-rmax,rmax))
    ax.set_zlim((-rmax,rmax))
    pl.draw()

    if eval_r(phi, a, e) < x_min:
        break

pl.show()

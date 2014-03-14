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
import mpl_toolkits.mplot3d as mpl3d

####################################
# Constants and Parameters

# general constants (CGS)
ggrav = 6.67e-8
c = 3.0e10
msun = 1.99e33
rsun = 6.96e10
direcs = [0, 1, 2] # indices corresponding to the three spatial directions

# simulation parameters
npoints = 1e4

# initial data
hulse_taylor_periastron = 1.1 * rsun
hulse_taylor_apastron = 4.8 * rsun
hulse_taylor_a = 0.5 * (hulse_taylor_periastron + hulse_taylor_apastron)
hulse_taylor_e = 1 - 2 / (hulse_taylor_apastron / hulse_taylor_periastron + 1)

# system parameters (CGS)
t = 0.0 # initial time
t_final = 1e5 # final time
init_major_axis = hulse_taylor_a
init_eccentricity = hulse_taylor_e
m1 = 1.4*msun # mass of body 1 (1.4 for typical neutron star)
m2 = 1.4*msun # mass of body 2
min_separation = 1e7 # simulation terminates if separation falls beneath this


####################################
# Helper Functions

def kron(i, j):
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
    Calculates the first time derivative of `quantity` at
    iteration `it`. If not enough data to compute, returns 0.
    """
    try:
        return (quantity[it] - quantity[it-1]) / dt
    except IndexError:
        return 0

def second_time_der(quantity, it):
    """
    Calculates the second time derivative of `quantity` at
    iteration `it`.  If not enough data to compute, returns 0.
    """
    try:
        return (quantity[it] - 2*quantity[it-1] + quantity[it-2]) / dt**2
    except IndexError:
        return 0

def third_time_der(quantity, it):
    """
    Calculates the third time derivative of `quantity` at
    iteration `it`.  If not enough data to compute, returns 0.
    """
    try:
        return (quantity[it] - 3*quantity[it-1] + 3*quantity[it-2]
                                                  - quantity[it-3]) / dt**3
    except IndexError:
        return 0

def eval_x(phi, a, e):
    """
    Calculates new binary separation vector.
    """
    return eval_r(phi, a, e) * np.array([np.cos(phi), np.sin(phi), 0.])

def eval_quad(x):
    """
    Calculates new reduced quadrupole tensor.
    """
    I_bar = np.zeros((3, 3))
    xnormsq = np.sum(x**2)
    for i in direcs:
        for j in direcs:
            I_bar[i,j] = mu * (x[i] * x[j] - kron(i,j) * xnormsq / 3.)
    return I_bar


def integrate_RK4(phi, a, e):
    """
    Integrates the true anomaly phi, the semi-major axis a, and the
    eccentricity e with the 4th-order Runge-Kutta method.
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
    """
    Calculates the binary separation distance r.
    From Equation 6 of the project description.
    """
    return a * (1 - e**2) / (1 + e * np.cos(phi))

def phi_RHS(r, a, e):
    """
    Righthand-side (time-derivative) of the true anomaly phi.
    From Equations 7 and 8 of the project description.
    """

    return np.sqrt(ggrav * M * a * (1 - e**2)) / r**2

def a_RHS(a, e):
    """
    Righthand-side (time-derivative) of semi-major axis a.
    From Equation 5.6 of Peters, 1964.
    """
    return ((-64./5.) * (ggrav**3 * m1 * m2 * M)
            * (1 + (73./24.) * e**2 + (37./96.) * e**4) 
            / (c**5 * a**3 * (1 - e**2)**3.5))

def e_RHS(a, e):
    """
    Righthand-side (time-derivative) of eccentricity e.
    From Equation 5.7 of Peters, 1964.
    """
    return ( (-304./15.) * e * (ggrav**3 * m1 * m2 * M)
             * (1 + (121./304.) * e**2) / (c**5 * a**4 * (1 - e**2)**2.5) )

def wave():
    """
    Function to calculate gravitational wave
    """
    pass

####################################
# main code body

pl.ion()

# set up initial conditions

rmax = 1.1 * init_major_axis * (1 + init_eccentricity) # 1.1 * apastron

M = m1 + m2
mu = (m1 * m2) / M

times = np.linspace(t, t_final, num=npoints)
dt = times[1] - times[0]

# array of true anomaly at each time step
phis = np.zeros(npoints)

# array of semi-major axis at each time step
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


# main loop 

for it, t in enumerate(times[:-1]):

    # Update quantities
    
    phis[it+1], axs[it+1], eccs[it+1] = integrate_RK4(phis[it], axs[it], eccs[it])
    xs[:,it+1] = eval_x(phis[it+1], axs[it+1], eccs[it+1])
    quads[:,:,it+1] = eval_quad(xs[:,it+1])
    
    r = np.linalg.norm(xs[:,it+1])
    hpluses[it+1] = second_time_der(quads[0,0,:] - quads[1,1,:], it+1) / r
    hcrosses[it+1] = 2 * second_time_der(quads[0,1,:], it+1) / r
    
    pos = np.zeros((2, 3))
    pos[0,:] =  (m2 / M) * xs[:,it+1]
    pos[1,:] = -(m1 / M) * xs[:,it+1]

    # Plot stuff
    
    pl.clf()
    fig = pl.gcf()
    ax = mpl3d.Axes3D(fig)
    ax.scatter(pos[:,0], pos[:,1], pos[:,2])
    ax.set_xlim((-rmax,rmax))
    ax.set_ylim((-rmax,rmax))
    ax.set_zlim((-rmax,rmax))
    pl.draw()

    if eval_r(phis[it+1], axs[it+1], eccs[it+1]) < min_separation:
        break

pl.show()

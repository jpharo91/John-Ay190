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
import cPickle as pickle
import time
import subprocess

####################################
# Constants and Parameters

# general constants (CGS units)
ggrav = 6.67e-8
c = 3.0e10
msun = 1.99e33
rsun = 6.96e10
direcs = [0, 1, 2] # indices corresponding to the three spatial directions

# initial data for Hulse-Taylor binary
hulse_taylor_periastron = 1.1 * rsun
hulse_taylor_apastron = 4.8 * rsun
hulse_taylor_a = 0.5 * (hulse_taylor_periastron + hulse_taylor_apastron)
hulse_taylor_e = 1 - 2 / (hulse_taylor_apastron / hulse_taylor_periastron + 1)

# parameters (CGS units)
movie = False # do we show a realtime movie of the system? much faster if False
movie_every = 10
output_every = 1000
timestamp_every = 10000
start_iteration = 19999000+1 # needs to be 0 if `input_dir` is `None`
input_dir = "inspiral2" # if None, a new simulation will be created
output_dir = "inspiral3" # name of directory to place output data in
kill_back_rxn = False # if False a and e will not be evolved

t_simulate = 8.5e-05 # total time to simulate
n_points = 1e7 # should be at least one point per 1000 s of simulated time
               # for phi evolution to be valid
init_a = 0.
init_e = 0.
m1 = 1.4*msun # mass of body 1 (1.4 for typical neutron star)
m2 = 1.4*msun # mass of body 2

if input_dir == None and start_iteration != 0:
    raise ValueError("If input_dir is None, then start_iteration must be 0")

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

def second_time_der(quantity0, quantity1, quantity2):
    """
    Calculates the second time derivative of a quantity based on values
    of that quantity at three times. The result is a
    forward derivative at the time of `quantity0`, a central derivative
    at the time of `quantity1`, and a backward derivative at the time of
    `quantity2`. To summarize the arguments:
    quantity0 = quantity(t_0)
    quantity1 = quantity(t_0 + dt)
    quantity2 = quantity(t_0 + 2*dt)
    (for some t_0)
    """
    return (quantity2 - 2*quantity1 + quantity0) / dt**2

def eval_r(phi, a, e):
    """
    Calculates the binary separation distance r.
    From Equation 6 of the project description.
    """
    return a * (1 - e**2) / (1 + e * np.cos(phi))

def eval_x(phi, a, e):
    """
    Calculates new binary separation vector.
    """
    return eval_r(phi, a, e) * np.array([np.cos(phi), np.sin(phi), 0.])

def eval_I_bar(x):
    """
    Calculates new reduced quadrupole tensor.
    """
    I_bar = np.zeros((3, 3))
    xnormsq = np.sum(x**2)
    for i in direcs:
        for j in direcs:
            I_bar[i,j] = mu * (x[i] * x[j] - kron(i,j) * xnormsq / 3.)
    return I_bar

def eval_h_plus(I_bar0, I_bar1, I_bar2, r):
    """
    Calculates plus-polarized gravitational-wave strain based on recent values
    of the reduced quadrupole I_bar.
    """
    return second_time_der(I_bar0[0,0] - I_bar0[1,1],
                           I_bar1[0,0] - I_bar1[1,1],
                           I_bar2[0,0] - I_bar2[1,1]) / r

def eval_h_cross(I_bar0, I_bar1, I_bar2, r):
    """
    Calculates cross-polarized gravitational-wave strain based on recent values
    of the reduced quadrupole I_bar.
    """
    return 2 * second_time_der(I_bar0[0,1], I_bar1[0,1], I_bar2[0,1]) / r

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

def RHS(quants):
    """
    Provides a length-3 array of the righthand-sides (time-derivatives)
    of phi, a, and e for use in integration schemes. `quants` is a length-3
    sequence with first element phi, second element a, and third element e.
    """
    phi, a, e = quants[0], quants[1], quants[2]
    r = eval_r(phi, a, e)
    return np.array([phi_RHS(r, a, e), a_RHS(a, e), e_RHS(a, e)])

def integrate_RK4(phi, a, e):
    """
    Integrates the true anomaly `phi`, the semi-major axis `a`, and the
    eccentricity `e` over one time step with the 4th-order Runge-Kutta method. 
    If the kill_back_rxn parameter is set to True, then `a` and `e` are not
    integrated.
    """
    old = np.array([phi, a, e])
    k1 = RHS(old)
    k2 = RHS(old + 0.5 * k1)
    k3 = RHS(old + 0.5 * k2)
    k4 = RHS(old + k3)
    new = old + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.
    if kill_back_rxn:
        new[1:] = old[1:] # set a and e to their old values
    return (new[0], new[1], new[2])

####################################
# main code body

start_time = time.time()

pl.ion()

if output_dir != None:

    subprocess.call("mkdir -p " + output_dir, shell=True)

    f = open("{}/parameters.txt".format(output_dir), 'w')
    f.write("t_simulate = {}\nn_points = {}".format(t_simulate, n_points))
    f.close()

if input_dir == None:
    # set up initial conditions (not continuing from a checkpoint)

    M = m1 + m2 # total mass
    mu = (m1 * m2) / M # reduced mass
    t = 0 # the current time
    phi = 0 # the true anomaly (orbital phase angle)
    a = init_a # the semi-major axis
    e = init_e # the eccentricity
    x = eval_x(phi, a, e) # the binary separation vector
    # I_bar2 is the current quadrupole, I_bar1 is what it was one time-step 
    # ago, and I_bar0 is what it was two time-steps ago
    I_bar0 = I_bar1 = I_bar2 = eval_I_bar(x)
    h_plus = 0 # the plus-polarized gravitational-wave strain
    h_cross = 0 # the cross-polarized gravitational-wave strain


    if output_dir != None:
        # initial data output
        data = {'t':t, 'phi':phi, 'a':a, 'e':e, 'x':x,
                'I_bar0':I_bar0, 'I_bar1':I_bar1, 'I_bar2':I_bar2,
                'h_plus':h_plus, 'h_cross':h_cross, 'm1':m1, 'm2':m2}
        f = open("{}/iteration-{}.pickle".format(output_dir, 0), 'wb')
        pickle.dump(data, f)
        f.close()


else:

    f = open("{}/iteration-{}.pickle".format(input_dir,start_iteration-1),'rb')
    data = pickle.load(f)
    f.close()
    
    t, phi, a, e, x = data['t'], data['phi'], data['a'], data['e'], data['x']
    h_plus, h_cross = data['h_plus'], data['h_cross']
    try:
        I_bar0, I_bar1, I_bar2 = data['I_bar0'], data['I_bar1'], data['I_bar2']
    except KeyError: # backwards compatibility with old format of output data
        I_bar0 = I_bar1 = I_bar2 = data['I_bar']
    try:
        m1, m2 = data['m1'], data['m2']
    except KeyError:
        pass # use the parameter values defined in this file
    M = m1 + m2 # total mass
    mu = (m1 * m2) / M # reduced mass

dt = t_simulate / n_points # the time-step
rmax = 1.1 * a * (1 + e) # 1.1 * apastron
final_iteration = start_iteration + n_points

# main loop 

for it in range(int(start_iteration), int(final_iteration)):

    # Time stamp

    if it % timestamp_every == 0:
        elapsed = time.time() - start_time
        print "Iteration {}/{:.0f} -- {}% -- t: {} -- a: {} -- Time elapsed: {:.0f} seconds ({:.1f} minutes) ({:.2f} hours)".format(
                                                    it, final_iteration,
                                                    (it-start_iteration)*100. / n_points,
                                                    t, a, elapsed, elapsed/60,
                                                            elapsed/3600)
 
    # Update quantities
    
    t += dt

    phi, a, e = integrate_RK4(phi, a, e)
    phi %= 2 * np.pi # keeps phi in [0, 2*pi)

    x = eval_x(phi, a, e)
    r = np.linalg.norm(x)

    I_bar0 = I_bar1
    I_bar1 = I_bar2
    I_bar2 = eval_I_bar(x)

    h_plus = eval_h_plus(I_bar0, I_bar1, I_bar2, r)
    h_cross = eval_h_cross(I_bar0, I_bar1, I_bar2, r)


    if movie and it % movie_every == 0:
    # Plot positions of binary stars

        pos = np.zeros((2, 3))
        pos[0,:] =  (m2 / M) * x
        pos[1,:] = -(m1 / M) * x
    
        pl.clf()
        fig = pl.gcf()
        ax = mpl3d.Axes3D(fig)
        ax.scatter(pos[:,0], pos[:,1], pos[:,2])
        ax.set_xlim((-rmax,rmax))
        ax.set_ylim((-rmax,rmax))
        ax.set_zlim((-rmax,rmax))
        pl.draw()


    if output_dir != None and it % output_every == 0:

        # Output data
        data = {'t':t, 'phi':phi, 'a':a, 'e':e, 'x':x,
                'I_bar0':I_bar0, 'I_bar1':I_bar1, 'I_bar2':I_bar2,
                'h_plus':h_plus, 'h_cross':h_cross, 'm1':m1, 'm2':m2}
        f = open("{}/iteration-{}.pickle".format(output_dir, it), 'wb')
        pickle.dump(data, f)
        f.close()

if movie:
    pl.show()

print "Done\a\a\a\a\a" # make noise when we're done

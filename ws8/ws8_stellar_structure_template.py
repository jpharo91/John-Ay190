#!/usr/bin/env python

import numpy as np
import scipy as sp


# global constants
ggrav = 6.67e-8
msun  = 1.99e33

# EOS parameters
# for white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG


#######################################
# function definitions
def tov_RHS(rad,rho,m):
    
    # RHS function
    
    rhs = np.zeros(2)
    if(rad > 1.0e-10):
        rhs[0] = (-1.0)*ggrav*m*rho/(rad*rad)
        rhs[1] = 4*np.pi*rho*rad*rad
    else:
        rhs[0] = 0.0
        rhs[1] = 0.0

    return rhs

def P2rho(pressure):
    return (pressure/polyK)**(1.0/polyG)

def tov_integrate_FE(rad,dr,p,rho,m):

    # Forward-Euler Integrator

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m

    # forward Euler integrator
    new = old + dr*tov_RHS(rad,p,rho,m)
    
    # assign outputs
    pnew = new[0]
    mnew = new[1]
    
    return (pnew,mnew)

def tov_integrate_RK2(rad, dr, p, rho, m):
    # Runge-Kutta 2

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m
    k2 = np.zeros(2)

    # Runge-Kutta 2 Integrator
    k1 = dr*tov_RHS(rad,rho,m)
    k2[0] = dr * tov_RHS(rad+.5*dr,P2rho(p+.5*k1[0]),m+.5*k1[0])[0]
    k2[1] = dr * tov_RHS(rad+.5*dr,P2rho(p+.5*k1[0]),m+.5*k1[1])[1]
    new = old + k2

    #assign outputs
    pnew = new[0]
    mnew = new[1]

    return (pnew, mnew)

def tov_integrate_RK3(rad, dr, p, rho, m):
    # Runge-Kutta 3

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m
    k2 = np.zeros(2)
    k3 = np.zeros(2)

    # Runge-Kutta 3 Integrator
    k1 = dr*tov_RHS(rad,rho,m)
    k2[0] = dr * tov_RHS(rad+.5*dr,P2rho(p+.5*k1[0]),m+.5*k1[1])[0]
    k2[1] = dr * tov_RHS(rad+.5*dr,P2rho(p+.5*k1[0]),m+.5*k1[1])[1]
    k3[0] = dr * tov_RHS(rad+dr,P2rho(p-k1[0]+2.0*k2[0]),m-k1[0]+2.0*k1[1])[0]
    k3[1] = dr * tov_RHS(rad+dr,P2rho(p-k1[0]+2.0*k2[0]),m-k1[0]+2.0*k1[1])[1]
    new = old + (1.0/6.0)*(k1 + 4.0*k2 + k3)

    #assign outputs
    pnew = new[0]
    mnew = new[1]

    return (pnew, mnew)

def tov_integrate_RK4(rad, dr, p, rho, m):
    # Runge-Kutta 4

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m
    k2 = np.zeros(2)
    k3 = np.zeros(2)
    k4 = np.zeros(2)

    # Runge-Kutta 4 Integrator
    k1 = dr*tov_RHS(rad,rho,m)
    k2[0] = dr * tov_RHS(rad+.5*dr,P2rho(p+.5*k1[0]),m+.5*k1[1])[0]
    k2[1] = dr * tov_RHS(rad+.5*dr,P2rho(p+.5*k1[0]),m+.5*k1[1])[1]
    k3[0] = dr * tov_RHS(rad+dr,P2rho(p-k1[0]+2.0*k2[0]),m-k1[0]+2.0*k1[1])[0]
    k3[1] = dr * tov_RHS(rad+dr,P2rho(p-k1[0]+2.0*k2[0]),m-k1[0]+2.0*k1[1])[1]
    k4[0] = dr * tov_RHS(rad+dr,P2rho(p+k3[0]),m)[0]
    k4[1] = dr * tov_RHS(rad+dr,P2rho(p+k3[0]),m+k3[1])[1]
    new = old + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)

    #assign outputs
    pnew = new[0]
    mnew = new[1]

    return (pnew, mnew)

#######################################

# set up grid
npoints = 1000
radmax = 2.0e8 # 2000 km
radius = np.linspace(0, radmax, num=npoints)
dr = radius[1]-radius[0]

# set up variables
press = np.zeros(npoints)
rho   = np.zeros(npoints)
mass  = np.zeros(npoints)

# set up central values
rho[0]   = 1.0e10
press[0] = polyK * rho[0]**polyG
mass[0]  = 0.0

# set up termination criterion
press_min = 1.0e-10 * press[0]

nsurf = 0
for n in range(npoints-1):
    
    (press[n+1],mass[n+1]) = tov_integrate_RK4(radius[n],
                                              dr,
                                              press[n],
                                              rho[n],mass[n])
    # check for termination criterion
    if(press[n+1] < press_min and nsurf==0):
        nsurf = n

    if(n+1 > nsurf and nsurf > 0):
        press[n+1] = press[nsurf]
        rho[n+1]   = rho[nsurf]
        mass[n+1]  = mass[nsurf]

    # invert the EOS to get density
    rho[n+1] = (press[n+1]/polyK)**(1.0/polyG)


print radius[nsurf]/1.0e5
print mass[nsurf]/msun




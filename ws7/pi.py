#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import numpy.random as random

# sample starting values
radius = 5.0
N = 1000000

# Function to calculate pi via Monte Carlo
def pi(radius, N):
    rad = radius * radius
    x = 2.0 * radius #Calculate sides of box in which circle in inscribed
    y = x
    n = 0 #Set initial counts to 0

    # Generate points in box
    xs = random.uniform(-x/2.0, x/2.0, N)
    ys = random.uniform(-y/2.0, y/2.0, N)

    # Check if points outside of box
    for i in range(len(xs)):
        xi = xs[i]
        yi = ys[i]
        x2y2 = xi*xi + yi*yi
        if x2y2 <= rad:
            n += 1

    pi = (4.0*n)/N
    return pi

# Function to test convergence at different values of N
def test(radius, Ns):
    results = []
    for N in Ns:
        val = pi(radius, N)
        frac = np.pi/val
        results.append((val, frac, N))
    return results

Ns = [10,100,1000,10000,100000,1000000]

convergence = test(radius, Ns)
print convergence

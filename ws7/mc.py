#!/usr/bin/env python

import numpy as np
import numpy.random as random

# The Function to be integrated
def f(x):
    return x*x + 1

# Function to integrate a function over a given interval via MC
def integrate(f, a, b, N):
    fb = f(b)
    fa = f(a)
    A1 = (b - a)*fb      # Calculate large box
    A2 = (b-a)*(fb-fa)   # Calculate box above f(a)
    xs = random.uniform(a, b, N) #Generate random points in A2
    ys = random.uniform(fa, fb, N)
    n = 0

    # Check if point less than function
    for i in range(len(xs)):
        xi = xs[i]
        yi = ys[i]
        if yi <= f(xi):
            n += 1

    print n
    value = A2*(n/np.float32(N)) + (A1 - A2)
    return value

# Function to check for convergence at various N
def test(a, b, Ns):
    results = []
    analytic = 22.0/3.0
    for N in Ns:
        val = integrate(f, a, b, N)
        frac = analytic/val
        results.append((val, frac, N))
    return results

Ns = [10,100,1000,10000,100000,1000000]

convergence = test(2.0, 3.0, Ns)
print convergence

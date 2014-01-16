#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as pl

# Helper function defining f(x)
def f(x):
    return (x*x*x - 5 * x*x + x)

# Helper function defining f'(x)
def f_prime(x):
    return (3 * x*x - 10 * x + 1)

# Function to calculate the forward difference approximation
def for_dif(func, begin, end, h):
    disc_range = []
    for i in np.arange(begin, end, h):
        disc_range.append((func(i + h) - func(i))/h)
    return disc_range

# Function to calculate the central difference approximation
def cen_dif(func, begin, end, h):
    disc_range = []
    for i in np.arange(begin, end, h):
        disc_range.append((func(i + h) - func(i - h))/(2*h))
    return disc_range

# Select resolution 
h1 = 1.0
h2 = h1/2.0

# Generate x values
disc_range = np.arange(-2, 6.2, h1)

derivatives = f_prime(disc_range) # Calculate analytic derivatives

# Calculate differences and absolute errors for both resolutions
forward_first = for_dif(f, -2, 6.2, h1)
central_first = cen_dif(f, -2, 6.2, h1)

y1 = forward_first - derivatives
y2 = central_first - derivatives

disc_range2 = np.arange(-2, 6.1, h2)
derivatives2 = f_prime(disc_range2)

forward_second = for_dif(f, -2, 6.1, h2)
central_second = cen_dif(f, -2, 6.1, h2)

y3 = forward_second - derivatives2
y4 = central_second - derivatives2

p1, = pl.plot(disc_range, y1, "r", linewidth=2)
p2, = pl.plot(disc_range, y2, "b", linewidth=2)
p3, = pl.plot(disc_range2, y3, "g", linewidth=2)
p4, = pl.plot(disc_range2, y4, "m", linewidth=2)

pl.xlim(min(disc_range), max(disc_range[:-1]))
pl.ylim(min(y1*1.05), max(y1*1.05))

pl.xlabel("X")
pl.ylabel("Difference")

pl.legend( (p1,p2,p3,p4), ("y1", "y2", "y3", "y4"), loc = 2, frameon=False  )

#pl.show()

pl.savefig('FinDifApprox.pdf')
    


        

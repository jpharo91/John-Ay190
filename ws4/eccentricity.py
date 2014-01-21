#!/usr/bin/env python

import numpy as np
import math as m

# Useful Parameters
T = 365.25635 #period in days
a = 1.496 * 10**6 #semi-major axis in km
e = 0.0167 #eccentricity
b = m.sqrt((1 - e*e)*(a*a))

# Helper function of Equation 4
def eq4(E):
    return E - (2.0*m.pi/T)*t - e*np.sin(E)

# Helper function of the derivative of Equation 4
def deriv(E):
    return 1.0 - e*np.cos(E)

# Function for computing Newton's method
def newton(func, prime, init, threshold = 0.1, max_iters = 10000):
    diff = float("inf")
    count = 0
    xi = init
    while (np.abs(diff) > threshold and count < max_iters):
        diff = float(func(xi)/prime(xi))
        xi -= diff
        count += 1
    print "# Iterations: " + str(count)
    return xi

# Function to calculate x
def x(E):
    return a*np.cos(E)

# Function to calculate y
def y(E):
    return b*np.sin(E)

# I used these to check my algorithm
#def para(x):
#    return x*x

#def para_prime(x):
#    return 2.0*x

# Let's take an initial guess of x=0
# with a threshold of 10^-10 to get the needed fractional accuracy

t = 91

print "For time t = " + str(t) + ", we have"
root = newton(eq4, deriv, 0, 10**(-10))
print "E = " + str(root)
X = x(root)
Y = y(root)
print "X = " + str(X)
print "Y = " + str(Y)

t = 182

print "For time t = " + str(t) + ", we have"
root = newton(eq4, deriv, 0, 10**(-10))
print "E = " + str(root)
X = x(root)
Y = y(root)
print "X = " + str(X)
print "Y = " + str(Y)

t = 273

print "For time t = " + str(t) + ", we have"
root = newton(eq4, deriv, 0, 10**(-10))
print "E = " + str(root)
X = x(root)
Y = y(root)
print "X = " + str(X)
print "Y = " + str(Y)

# Now suppose the eccentricity changes
e = 0.99999
t = 91

print "For time t = " + str(t) + ", we have"
root = newton(eq4, deriv, 0, 10**(-10))
print "E = " + str(root)
X = x(root)
Y = y(root)
print "X = " + str(X)
print "Y = " + str(Y)

# Let's try making a better guess informed by our knowledge of the geometry 
# of the situation

print "\nTry again with an initial guess of pi"
print "For time t = " + str(t) + ", we have"
root = newton(eq4, deriv, np.pi, 10**(-10))
print "E = " + str(root)
X = x(root)
Y = y(root)
print "X = " + str(X)
print "Y = " + str(Y)

t = 182

print "For time t = " + str(t) + ", we have"
root = newton(eq4, deriv, np.pi, 10**(-10))
print "E = " + str(root)
X = x(root)
Y = y(root)
print "X = " + str(X)
print "Y = " + str(Y)

t = 273

print "For time t = " + str(t) + ", we have"
root = newton(eq4, deriv, np.pi, 10**(-10))
print "E = " + str(root)
X = x(root)
Y = y(root)
print "X = " + str(X)
print "Y = " + str(Y)



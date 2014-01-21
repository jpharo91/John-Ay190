#!/usr/bin/env python

import numpy as np
import math as m

def midpoint(func, a, b, step):
    points = np.arange(a+step/2, b, step)
    mid_vals = func(points)
    return step * sum(mid_vals)

def trapezoid(func, a, b, step):
    points = np.arange(a, b+step, step)
    ends = func(points)
    return (step/2) * (2*sum(ends) - ends[0] - ends[-1])

def simpson(func, a, b, step):
    mids = np.arange(a+step/2, b, step)
    ends = np.arange(a, b+step, step)
    f_mids = func(mids)
    f_ends = func(ends)
    return (step/6)*(2*sum(f_ends) + 4*sum(f_mids) - f_ends[0] - f_ends[-1])

def error(analytic, numerical):
    absErr = analytic - numerical
    relErr = (analytic - numerical)/analytic
    return (absErr, relErr)

# The analytic solution for the integral of sin(x) from 0 to pi is 2.0

mid = midpoint(np.sin, 0, m.pi, .01)
print "\nMidpoint Rule"
print "Val: " + str(mid)
print "Errors: " + str(error(2.0, mid))

trap = trapezoid(np.sin, 0, m.pi, .01)
print "\nTrapezoidal Rule"
print "Val: " + str(trap)
print "Errors: " + str(error(2.0, trap))

simp = simpson(np.sin, 0, m.pi, .01)
print "\nSimpsons Rule"
print "Val: " + str(simp)
print "Errors: " + str(error(2.0, simp))

print "\nNow let's reduce the step size by a factor of 10"

mid = midpoint(np.sin, 0, m.pi, .001)
print "\nMidpoint Rule"
print "Val: " + str(mid)
print "Errors: " + str(error(2.0, mid))

trap = trapezoid(np.sin, 0, m.pi, .001)
print "\nTrapezoidal Rule"
print "Val: " + str(trap)
print "Errors: " + str(error(2.0, trap))

simp = simpson(np.sin, 0, m.pi, .001)
print "\nSimpsons Rule"
print "Val: " + str(simp)
print "Errors: " + str(error(2.0, simp))

# Now let's integrate the function g(x) = x*sin(x) on the same interval
# which has an analytic solution of pi

def g(x):
    return x*np.sin(x)

print "\nNew function"

mid = midpoint(g, 0, m.pi, .01)
print "\nMidpoint Rule"
print "Val: " + str(mid)
print "Errors: " + str(error(m.pi, mid))

trap = trapezoid(g, 0, m.pi, .01)
print "\nTrapezoidal Rule"
print "Val: " + str(trap)
print "Errors: " + str(error(m.pi, trap))

simp = simpson(g, 0, m.pi, .01)
print "\nSimpsons Rule"
print "Val: " + str(simp)
print "Errors: " + str(error(m.pi, simp))

print "\nNow let's reduce the step size by a factor of 10"

mid = midpoint(g, 0, m.pi, .001)
print "\nMidpoint Rule"
print "Val: " + str(mid)
print "Errors: " + str(error(m.pi, mid))

trap = trapezoid(g, 0, m.pi, .001)
print "\nTrapezoidal Rule"
print "Val: " + str(trap)
print "Errors: " + str(error(m.pi, trap))

simp = simpson(g, 0, m.pi, .001)
print "\nSimpsons Rule"
print "Val: " + str(simp)
print "Errors: " + str(error(m.pi, simp))



#!/usr/bin/env python

from numpy import float32

# establish constants from Equation 1
oneThird = float32(1.0/3)
thirteenThirds = float32(13.0/3)
fourThirds = float32(4.0/3)

# function to recursively calculate x_n as in Equation 1
def thirds(n):
    terms = [float32(1),oneThird]
    for i in range (2,n):
        terms.append(thirteenThirds*terms[-1] - fourThirds * terms[-2])
    return terms

# calculate actual values for error comparison
def precise(n):
    return [(1.0/3) ** x for x in range(n)] # get maximum precision

# function to calculate the absolute error
def absError(n):
    absErr = []
    prec = precise(n)
    recursive = thirds(n)
    for i in range(n):
        absErr.append(recursive[i] - prec[i])
    return absErr

# function to calculate the relative error
def relError(n):
    relErr = []
    prec = precise(n)
    recursive = thirds(n)
    for i in range(n):
        relErr.append((recursive[i] - prec[i])/prec[i])
    return relErr

# function to print results in a format useful for copying into LaTeX
def latex(n, vals, absErr, relErr):
    final = ''
    for i in range(n):
        row = [i, vals[i], absErr[i], relErr[i]]
        for j, item in enumerate(row):
            row[j] = str(item)
        row = ' & ' .join(row) + ' \\\\ \n'
        final += row
    return final

values = thirds(16)
absolute = absError(16)
relative = relError(16)

results = latex(16, values, absolute, relative)
print results



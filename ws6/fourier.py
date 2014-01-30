#!/usr/bin/env python

import numpy as np
from timeit import timeit
import matplotlib.pyplot as pl

# The square root of negative one is i, not j
i = 1j

# Function to compute a discrete Fourier Transform
def dft(x):
    N = len(x)
    w = np.exp(-(2*np.pi*i)/N)
    theMatrix = []
    for j in range(N):
        row = []
        for k in range (N):
            row.append(w ** (j*k))
        theMatrix.append(row)
    return np.dot(theMatrix, x)

# Abandoned attempt at faster version using map
#def dft2(x):
#    N = len(x)
#    w = np.exp(-(2*np.pi*i)/N)
#    theMatrix = []
#    for j in range(N):
#        k
#        row = map(w **)

x = np.random.randn(10)

# Verify that my function gets the same results as numpy's
print dft(x)
print np.fft.fft(x)
# it does

xs = range(10,101)
dft_times = []
fft_times = []

# Time the two dft functions
for i in xs:
    dft_times.append(timeit("dft(x)", number=10, \
        setup="from fourier import dft; import pylab; x=pylab.randn(%d)" % i))
    fft_times.append(timeit("np.fft.fft(x)", number=10, \
        setup="import numpy as np; import pylab; x=pylab.randn(%d)" % i))

p1, = pl.semilogy(xs, dft_times, "r", linewidth=2)
p2, = pl.semilogy(xs, fft_times, "b", linewidth=2)

pl.xlabel("N")
pl.ylabel("Log Time (s)")

pl.legend( (p1,p2), ("My Function","Numpy's Function"), loc=2, frameon=False )

#pl.show()

pl.savefig("Times.pdf")

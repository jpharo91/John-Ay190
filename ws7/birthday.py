#!/usr/bin/env python

import numpy.random as random
import numpy as np
import matplotlib.pyplot as pl

# Function to check for matching birthdays in group of given size
def matches(size, N):
    birthdays = []
    n = 0

    for test in range(N):
        taken = {} #Place to put taken dates

        for s in range(size):
            bday = random.randint(0, 365)
            if bday in taken:
                n += 1
                break
            taken[bday] = 1 #Add to taken dates

    frac = n/np.float32(N)
    return frac

# Function to check for matches with N=10000
def birthdays(size):
    return matches(size, 10000)

# Function to look at convergence for different values of N
def convergence(size):
    Ns = [10,100,1000,10000,100000,1000000] #Wide range of N's
    results = []
    for N in Ns:
        val = matches(size, N)
        frac = .5/val
        results.append((val, frac, N))
    return results

sizes = range(2, 51)

probs = map(birthdays, sizes)
ys = [.5 for number in range(0, 51)]

p1, =pl.plot(sizes, probs, "bo", linewidth=2)
p2, =pl.plot(range(0, 51), ys, "r", linewidth=2)

pl.xlabel("Number of People")
pl.ylabel("Probability of a Shared Birthday")

pl.legend( (p1,p2), ("Probability","Minimum"), loc=2, frameon=False )

pl.savefig("Birthday.pdf")

#Now that we know the answer is 23, check convergence to P=0.5
test = convergence(23)
print test

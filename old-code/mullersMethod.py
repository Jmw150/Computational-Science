#!/usr/bin/env python

# Muller's method
# made by: Jordan Winkler
# finds an approximation of a root using Muller's method
# This program is weak against injection attacks.
# option -f for fractions, default floating point
#        -d for debugging or verbose mode
#        -i [number] for iterations of function
#        -p [number] for precision

# could have been defined recursively or with a do while loop
# f is a function, the rest are numbers
# returns number on success, string on failure
def mullersMethod (f, p0, p1, p2, iterations, tolerance) :
    # set up parabola interpolation
    h1 = p1 - p0
    h2 = p2 - p1
    s1 = (f(p1) - f(p0))/float(h1)
    s2 = (f(p2) - f(p1))/float(h2)
    d  = (s2 - s1)/float((h2 + h1))
    i = 3

    while i <= iterations :
        # calc base
        b = s2 + h2*d
        D = ((b**2 - 4*f(p2)*d)+0j)**0.5 #convert to complex 
        
        # pick the larger base
        if abs(b-D) < abs(b+D) :
            E = b + D
        else :
            E = b - D

        # calc rest of equation
        h = -2*f(p2)/E
        p = p2 + h
        
        # if within tolerance
        if abs(h) < tolerance :
            return p

        # update values to go again
        p0 = p1
        p1 = p2
        p2 = p
        h1 = p1 - p0
        h2 = p2 - p1
        s1 = (f(p1) - f(p1))/(h1)
        s2 = (f(p2) - f(p1))/(h2)
        d  = (s2 - s1)/(h2 + h1)
        i += 1
    
    return " error "
        


# libraries 
import sys
import sympy as sp
from numpy import *
from fractions import Fraction

fstring = sys.argv[1]
f = lambda x : eval(fstring)
p0 = eval(sys.argv[2])
p1 = eval(sys.argv[3])
p2 = eval(sys.argv[4])
iterat = 50
tol = 0.0001

print (mullersMethod(f, p0, p1, p2, iterat, tol))

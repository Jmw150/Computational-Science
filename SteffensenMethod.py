#!/usr/bin/env python

# finds an approximation of a root using Newton's Method
# default iterations is 15

from sympy import *
import sys

if len(sys.argv) != 3 :
    print("steffensen [function] [seed]")
    exit()

# function:string, seed:number
def steffensenMethod (function, seed, precision) :

    # error check
    if precision == 0 :
        exit()
    
    # setup
    h = precision
    x = seed
    f = lambda x : eval(function)
    g = lambda x : (f(x+h)-f(x))/h
    iterations = 50

    while (iterations > 0) :
        x -= f(x)/float(g(x))
        iterations -= 1
    return x
    
function = sys.argv[1]
seed = eval(sys.argv[2])

# solution close to
print ((steffensenMethod(function,seed,0.1)))

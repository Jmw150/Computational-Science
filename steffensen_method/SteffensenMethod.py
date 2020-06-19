#!/usr/bin/env python

# finds an approximation of a root using steffensen's Method
# default iterations is 50

import sys

if len(sys.argv) != 3 :
    print("steffensen [function] [seed]")
    exit()

# function:string, seed:number
def steffensenMethod (function, seed, precision, iterations) :

    # error check
    if precision == 0 :
        exit()
    
    # setup
    i = iterations
    h = precision
    x = seed
    f = lambda x : eval(function)
    g = lambda x : (f(x+h)-f(x))/h

    while (i > 0) :
        x -= f(x)/float(g(x))
        i -= 1
    return x
    
function = sys.argv[1]
seed = eval(sys.argv[2])

# solution close to
print ((steffensenMethod(function,seed,0.1,50)))

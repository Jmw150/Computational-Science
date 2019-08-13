#!/usr/bin/env python

# finds an approximation of a root using Newton's Method
# default iterations is 15

from sympy import *
import sys

if len(sys.argv) != 3 :
    print("newtonsMethod [function] [seed]")
    exit()

# function:string, seed:number
def newtonsMethod (function, seed) :
    x = Symbol('x')
    functionp = str(diff(function,x)) #sympy part

    x = seed
    f = lambda x : eval(function)
    fp = lambda x : eval(functionp)
    iterations = 50
    while (iterations > 0) :
        x -= f(x)/float(fp(x))
        iterations -= 1
    return x
    
function = sys.argv[1]
seed = eval(sys.argv[2])

# solution close to
print ((newtonsMethod(function,seed)))

#!/usr/bin/env python

# Newton's method
# made by: Jordan Winkler
# finds an approximation of a root using Newton's Method
# This program is weak against injection attacks.
# option -f for fractions, default floating point
#        -d for debugging or verbose mode
#        -i [number] for iterations of function
#        -p [number] for precision


# function:string, seed:number
def newtonsMethod (function, seed, iterations) :
    x = sp.Symbol('x')
    functionp = str(sp.diff(function,x)) #sympy part
    
    if argexists(fractionArg) :
        x = Fraction(seed)
    else :
        x = seed
    f = lambda x : eval(function)
    fp = lambda x : eval(functionp)
   
    if (argexists(fractionArg)) :
        zeroIsh = nearzero
    else :
        zeroIsh = 0
    
    # actual algorithm
    while (iterations > 0 and (f(x) < -zeroIsh or f(x) > zeroIsh)) :
        if (argexists(debug)) :
            print("x: " + str(x))
            print("i: " + str(iterations))
        if argexists(fractionArg) :
            x -= f(x)/(fp(x)) #can get crazier fractions each cycle
        else :
            x -= f(x)/float(fp(x))
        iterations -= 1
    return x

##################################################################
# argexists [string]
#
# checks if argument from command line was given
#
# impure function:
#     input: string
#     output: location if argument was given, else 0
##################################################################
def argexists (string) :
    length = len(sys.argv)
    i = 1
    while i < length :
        if sys.argv[i] == string :
            return i
        i += 1
    return 0

# command line options
debug = "-d"
fractionArg = "-f"
iteration = "-i"
precision = "-p"

# libraries 
import sys
import sympy as sp
from numpy import *
from fractions import Fraction

if len(sys.argv) == 1 :
    print("newtonsMethod [function] [seed]"
          "\n"
          "options:\n"
          "-f for fractions, default floating point\n"
          "-d for debugging or verbose mode\n"
          "-i [number] for iterations of function\n"
          "-p [number] for precision\n")
    exit()

# iterations
if argexists(iteration) :
    iterations = eval(sys.argv[(argexists(iteration)+1)])
else :
    iterations = 10 #converges quickly, fractions get hard to compute quickly too

# Precision of Bisection Method? (fractions only) 
if argexists(fractionArg) :
    if argexists(precision) :
        nearzero = Fraction(1,eval(sys.argv[(argexists(precision)+1)]))
    else :
        nearzero = Fraction(1,2**8) 
    
function = sys.argv[1]
if argexists(fractionArg) :
    seed = Fraction(eval(sys.argv[2]))
else :
    seed = eval(sys.argv[2])

# solution close to
print ((newtonsMethod(function, seed, iterations)))

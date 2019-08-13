#!/usr/bin/env python

# Steffensen's method
# made by: Jordan Winkler
# finds an approximation of a root using Steffesen's Method
# This program is weak against injection attacks.
# option -f for fractions, default floating point
#        -d for debugging or verbose mode
#        -i [number] for iterations of function
#        -p [number] for precision


# function:string, seed:numbers, iterations:numbers
def steffensensMethod (function, seed, iterations) :
    f = lambda x : eval(function)

    x = seed

    tolerance = nearzero

    while (iterations > 0 and (f(x) > tolerance or f(x) < -tolerance)) :
        if argexists(debug) :
            print ("x = " + str(x))
            print ("i = " + str(iterations))
        if argexists(fractionArg) :
            x -= f(x)/((f(x+f(x)) - f(x))/(f(x)))
        else :
            x -= float(f(x))/((f(x+f(x)) - f(x))/(f(x)))
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

#if len(sys.argv) != 3 :
#    print("newtonsMethod [function] [seed]")
#    exit()

# iterations
if argexists(iteration) :
    iterations = eval(sys.argv[(argexists(iteration)+1)])
else :
    iterations = 10 #converges quickly, fractions get hard to compute quickly too

# precision
if argexists(precision) :
    nearzero = Fraction(1,eval(sys.argv[(argexists(precision)+1)]))
else :
    if argexists(fractionArg) :
        nearzero = Fraction(1,2**8) #more precise fractions are huge
    else :
        nearzero = Fraction(1,2**32) 
    
function = sys.argv[1]
if argexists(fractionArg) :
    seed = Fraction(eval(sys.argv[2]))
else :
    seed = eval(sys.argv[2])

# solution close to
print ((steffensensMethod(function, seed, iterations)))

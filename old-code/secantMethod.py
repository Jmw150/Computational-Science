#!/usr/bin/env python

# secant Method
# made by: Jordan Winkler
# finds an approximation of a root using secant Method
# This program is weak against injection attacks.
# option -f for fractions, default floating point
#        -d for debugging or verbose mode
#        -i [number] for iterations of function
#        -p [number] for precision


# function:string, seed:numbers, iterations:numbers
def secantMethod (function, seed1, seed2, iterations) :
    f = lambda x : eval(function)

    # these probably need to be flipped
    x = seed1
    x_t = seed2 # x's tail

    tolerance = nearzero

    while (iterations > 0 and (abs(x-x_t) > tolerance)) :
        if argexists(debug) :
            print ("x       = " + str(x))
            print ("x_t     = " + str(x_t))
            print ("i       = " + str(iterations))
        temp = x
        if argexists(fractionArg) :
            x -= f(x)/((f(x) - f(x_t))/(x-x_t))
        else :
            x -= f(x)/float((f(x) - f(x_t))/(x-x_t))
        x_t = temp
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
    nearzero = Fraction(1,2**8) 
    
function = sys.argv[1]
if argexists(fractionArg) :
    seed1 = Fraction(eval(sys.argv[2]))
    seed2 = Fraction(eval(sys.argv[3]))
else :
    seed1 = eval(sys.argv[2])
    seed2 = eval(sys.argv[3])

# solution close to
print ((secantMethod(function, seed1, seed2, iterations)))

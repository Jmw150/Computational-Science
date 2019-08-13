#!/usr/bin/env python

# bisection method
# made by: Jordan Winkler
# This is a demonstration of the bijection method. 
# This program is weak against injection attacks.
# option -f for fractions, default floating point
#        -d for debugging or verbose mode
#        -i [number] for iterations of function
#        -p [number] for precision


#functions

##################################################################
# bisection [function] [start] [stop] [iterations]
# 
# checks if the bisection method would work on a function of a given 
# domain
# 
# impure function: 
#     takes in a function and its domain
#     returns True/False if the bisection method would work
#     behavior changes if -f, or -d are terminal arguments
##################################################################
def bisection (f, start, stop, iterations) :
    if argexists(debug) :
        print ("||function start||") 
        print("function: " + functionstring)
        print("start: " + str(start))
        print("stop: " + str(stop))
        print((f(start)*f(stop)) < 0)
        print("iterations left: " + str(iterations))

    if argexists(fractions) :
        mid = (start+stop)/Fraction(2)
    else :
        mid = float((start+stop))/2

    if iterations == 0 :
        return mid

    if argexists(fractions) :
        if f(start) > nearzero :
            if f(mid) > nearzero :
                return bisection(f, mid, stop, iterations-1)
            elif f(mid) < -nearzero :
                return bisection(f, start, mid, iterations-1)
            elif f(mid) >= -nearzero and f(mid) <= nearzero :
                return mid
        if f(start) < -nearzero :
            if f(mid) < -nearzero :
                return bisection(f, mid, stop, iterations-1)
            elif f(mid) > nearzero :
                return bisection(f, start, mid, iterations-1)
            elif f(mid) >= -nearzero and f(mid) <= nearzero :
                return mid
    else :
        if f(start) > 0 :
            if f(mid) > 0 :
                return bisection(f, mid, stop, iterations-1)
            elif f(mid) < 0 :
                return bisection(f, start, mid, iterations-1)
            elif f(mid) == 0 :
                return mid
        if f(start) < 0 :
            if f(mid) < 0 :
                return bisection(f, mid, stop, iterations-1)
            elif f(mid) > 0 :
                return bisection(f, start, mid, iterations-1)
            elif f(mid) == 0 :
                return mid

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

##################################################################
# bisectable [function] [start] [stop]
# 
# checks if the bisection method would work on a function of a given 
# domain
# 
# pure function: 
#     takes in a function and its domain
#     returns True/False if the bisection method would work
##################################################################
def bisectable (f, start, stop) :
    return (f(start)*f(stop)) < 0

# command line definitions
debug = "-d"
fractions = "-f"
iterate = "-i"
precis = "-p"


#libraries
import sys
from math import * 
# to save a little time fractions are not linked unless used
if argexists(fractions) :
    from fractions import Fraction

# catch common errors before starting
if ((len(sys.argv) < 4 or len(sys.argv) > 10) or
    (len(sys.argv) == 5 and not (argexists(fractions) or argexists(debug))) or 
    (len(sys.argv) == 6 and not (argexists(iterate) or argexists(precis))) or
    (len(sys.argv) == 6 and not (argexists(fractions) and argexists(debug))) or
    (len(sys.argv) == 10 and not (argexists(fractions) or argexists(debug) or argexists(iterate) or argexists(precis)))) :
    print ("bisectionMethod [function] [start] [stop] (options)\n"
           "\n"
           "options:\n"
           "-f fractions\n"
           "-d debug mode\n"
           "-i [num] iterations\n"
           "-p [num] precision, only if using fractions\n"
          )
    exit()


# user input
functionstring = sys.argv[1]
f = lambda x : eval(functionstring)
if argexists(fractions) :
    start = Fraction(eval(sys.argv[2]))
    stop = Fraction(eval(sys.argv[3]))
else :
    start = eval(sys.argv[2])
    stop = eval(sys.argv[3])

# How many iterations of the Bisection Method?
if argexists(iterate) :
    iterations = eval(sys.argv[(argexists(iterate)+1)])
else :
    iterations = 100

# Precision of Bisection Method? (fractions only) 
if argexists(fractions) :
    if argexists(precis) :
        nearzero = Fraction(1,eval(sys.argv[(argexists(precis)+1)]))
    else :
        nearzero = Fraction(1,2**128)

#after getting user input this is the implementation 
if argexists(debug) :
    print("||global||")
    print(functionstring)
    print(start)
    print(stop)
    print(iterations)

if not bisectable(f, start, stop) :
    print "Error: This function is not ideal for the bisection method"
else :
    bisect = bisection(f, start, stop, iterations)
    print(str(bisect))

   

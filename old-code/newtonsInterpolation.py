#!/usr/bin/env python

# Newton's interpolation
# made by: Jordan Winkler
# finds the interpolling polynomial of a function 
# option 
#        -fun function instead of f(x) list
#        -d   degree of interpol
#        -x   x table here
#        -f   f(x) table here
#        -fx  give back value here (default is giving back string)
#        -xin x as an interval, and n size


usedFun = "-fun"
useDegree = "-d"
xtable = "-x"
ftable = "-f"
valueOnly = "-fx"
xInterval = "-xin"

import sys
from math import *
from fractions import Fraction

if len(sys.argv) == 1 :
    print("newtonsInterpolation\n\n"
          "options:\n"
          "-fun function instead of f(x) list\n"
          "-d   degree of interpol\n"
          '-x   x table here\n'
          '-f   f(x) table here\n'
          '-fx  give back value here (default is giving back string)\n'
          '-xin x as an interval, and n size\n')
    exit(1)

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

# x and f should be global
# this makes default location 1 unless specified
if argexists(xtable) :
    x = eval(sys.argv[argexists(xtable)+1])

# -xin [start,end] intervals, assumes function argument was also called
if argexists(xInterval) :
    x = [] #clean above
    startNend = eval(sys.argv[argexists(xInterval)+1])
    start = startNend[0]
    end = startNend[1]
    size = eval(sys.argv[argexists(xInterval)+2])
    increment = (end-start)/float(size)
    while start <= end :
        x.append(start)
        start += increment

if argexists(ftable) :
    f = eval(sys.argv[argexists(ftable)+1])

# option to use a function instead of f(x) values
if argexists(usedFun) :
    f = []
    function = lambda x : eval(sys.argv[argexists(usedFun)+1])
    index = 0
    while index < len(x) :
        f.append(function(x[index]))
        index += 1

# option to name the degree of this polynomial
if argexists(useDegree) :
    degree = eval(sys.argv[argexists(useDegree)+1])
else :
    degree = len(x)

# divided difference notation
def F(i,j) :
   if (j == 0) : 
       return float(f[i])
   #elif (j == 1) : # base case
   #    return float(f[i] - f[i-1])/float(x[i] - x[i-1])
   else :
       return (F(i,j-1) - F(i-1,j-1))/float(x[i] - x[i-j])

# creates a string of the newtons polynomial for later, and quicker,
# use.
# P(x) = F(0,0) + sum_{i=1 to n} (F(i,i)*product_{j=0 to i-1} (x-x[j]))
def makeNewtonPoly (degree=len(x)) :
    if degree < len(x)-1 :
        n = degree
    else :
        n = len(x) - 1
    P = str(F(0,0))

    i = 1
    while i <= n :
        P = P + "+(" + str(F(i,i)) + ")"
        i += 1
        
        j = 0
        while j < (i-1) : #did one more time than asked
            P = P + "*(x-(" + str(x[j]) + "))"
            j += 1

    return P


P = makeNewtonPoly(degree)
newPoly = lambda x : eval(P)


if argexists(valueOnly) :
    xval = eval(sys.argv[argexists(valueOnly)+1])
    print(newPoly(xval))
else :
    print (P)


# test function for class
def test () :
    maximum = abs(function(-1)-newPoly(-1)) # function not defined
    i = 1
    while i <= 40 :
        xt = -1 + 2*(i-1)/float(40)
        x  = -1 + 2*i/float(40)
        if abs(function(x)-newPoly(x)) > abs(function(xt)-newPoly(xt)) :
            maximum = abs(function(x)-newPoly(x))
        i += 1
    return maximum
        
#print(test()) # used in class, really should be written as a option


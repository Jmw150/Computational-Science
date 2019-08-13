#!/usr/bin/env python

# gets the length of a string described by a function

from sympy import *
import sys

if len(sys.argv) != 4 :
    print("arcLength [function] [start] [end]")
    exit()

# function is a string, rest are numbers
# output is a sympy object
def arcLength (function, start, end) :
    x = Symbol('x')
    functionp = str(diff(function,x))
    lengthStr = str(integrate("(1 + " + functionp + "**2)**(1/2)", x))
    lengthInt = lambda x : eval(lengthStr)
    return (lengthInt(end) - lengthInt(start))

ffunction = sys.argv[1]
fstart = eval(sys.argv[2])
fend = eval(sys.argv[3])

arcLengthEx = arcLength(ffunction, fstart, fend)

print ("arc length: " + str(arcLengthEx))


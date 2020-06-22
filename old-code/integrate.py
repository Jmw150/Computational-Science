#!/usr/bin/env python

from math import * # same as using include and namespace
from sympy import *
import sys

if len(sys.argv) != 2 :
     print("use one argument")
     exit()

x = Symbol('x') 

function = sys.argv[1]

# only handles string input, passes back sympy object
integral = integrate(function,x) 
# remove pointless sympy object, use string, not needed
integral = str(integral)

print(integral + " + C")

#x = 3
#
#print("if x=3 ")
#print(eval(integral))
#
#print("integral again")
#
#x = Symbol('x') 
#secondIntegral = integrate(integral,x)
#
#print(secondIntegral)
#

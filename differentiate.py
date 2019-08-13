#!/usr/bin/env python

from math import * # same as using include and namespace
from sympy import *
import sys

if len(sys.argv) != 2 :
     print("use one argument")
     exit()

# identify what is a variable
x = Symbol('x') 

# function = eval(raw_input("function "))

function = sys.argv[1]

# only handles string input, passes back sympy object
derivative = diff(function,x) 
# remove pointless sympy object, use string, not needed
derivative = str(derivative)

print(derivative)

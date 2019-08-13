#!/usr/bin/env python

# may actually not give accurate answers
# constant was off by a lot, x**1 was opposite sign

from math import * # same as using include and namespace
from sympy import *
from mpmath import fac #factorial
import sys

x = Symbol('x')

if len(sys.argv) != 4 :
    print ("taylorPolynomial [function] [n-derivatives] [centered on]")
    exit()

function = sys.argv[1]
k = nDerivatives = eval(sys.argv[2])
xnot = eval(sys.argv[3])

derivative = []
derivative.append(function)
index = 0
while index <= nDerivatives :
    derivative.append(str(diff(derivative[index],x)))
    index += 1

derivativef = []
index = 0
while index <= nDerivatives :
    derivativef.append(lambda x: eval(derivative[index])) # <3
    index += 1

factorials = []
index = 0
while index <= nDerivatives :
    factorials.append(fac(index))
    index += 1

polynomial = []
index = 0
while index <= nDerivatives :
    polynomial.append((x-xnot)**index)
    index += 1

taylorPolynomial = []
index = 0
while index <= nDerivatives :
    taylorPolynomial.append(derivativef[index](xnot)/factorials[index] * polynomial[index])
    index += 1

print (taylorPolynomial)

taylorPolynomialf = 0
index = 0
while index <= nDerivatives :
    if index % 2 == 0 :
        taylorPolynomialf += taylorPolynomial[index]
    else :
        taylorPolynomialf -= taylorPolynomial[index]
    index += 1

print (taylorPolynomialf)

taylorPolynomialfunction = lambda x : eval(str(taylorPolynomialf))

test = eval(raw_input("test value: "))

functionf = lambda x : eval(function)

print("T(x) " + str(taylorPolynomialfunction(test)))
print("f(x) " + str(functionf(test)))



    



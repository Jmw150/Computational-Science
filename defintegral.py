#!/usr/bin/env python

# lump sum of sum methods (and definite integrals)
# made by: Jordan Winkler
# not a lot of error handling went into this program
#
# common methods:
###################
#Rectangular rule: $\int^b_a f(a)\, dx = f(a)(b-a)$
#  error : $f'(me)(b-a)^2/2$
#
#Trapezoid Rule: $\int^b_a f(x)\, dx =$$ \f{b-a}{2}(f(a) + f(b))$
#  error : $-\f{f''(me)(b-a)^3}{12}$
#composite : $\int^b_a f(x) dx \approx \f{b-a}{n*2}(f(a) + 2 \sum^{n-1}_{i=1} f(x_i) + f(b))$
#  error : $-\f{1}{12} f''(me)(b-a) h^2$
#
#Simpson's Rule: $\int^b_a f(x) \, dx =$$ \f{b-a}{6} (f(a) + 4f(\f{a+b}{2}) + f(b))$
#  error : $-\f{1}{90} (\f{b-a}{2})^5 f^{(4)}(me)$
#composite : $\int^b_a f(x) \, dx \approx \f{b-a}{3*n}(f(x_0) + 2 \sum^{n/2 - 1}_{j=1} f(x_{2j}) + 4\sum^{n/2}_{j=1} f(x_{2j-1}) + f(x_n))$
#  error : $-\f{h^4(b-a)}{180} f^{(4)}(me)$
#
#


def compositeTrap (f,a,b,n) :
    h = (b-a)/float(n)

    # get sum part
    summ = 0
    i = 1
    while (i <= n-1):
        xi = a + i*h
        summ += f(xi)
        i += 1

    return h/2*(f(a) + 2*summ + f(b))

def compositeSimson (f,a,b,n) :
    h = (b-a)/float(n)

    # get sum part
    sum1 = 0
    j = 1
    while (j <= n/2-1):
        index = 2*j
        sum1 += f(a + index*h)
        j += 1

    sum2 = 0
    j = 1
    while (j <= n/2):
        index = 2*j-1
        sum2 += f(a+index*h)
        j += 1

    return h/3*(f(a) + 2*sum1 + 4*sum2 + f(b))

# need to get n high enough to put error bound under user given bound
# turn h = (a-b)/n
# solve for n
#follow by hand algorithm from previous problems
   #error : $-\f{h^4(b-a)}{180} f^{(4)}(me)$ = error_bound
#function from a to b, and error desired
def comSimErr(f, a, b, e) :
    return 5 #'stub'

def comTrapErr(e) :
    return 5 #'stub'

import sys
from math import *
from fractions import Fraction

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

# function as a set of ordered pairs
# set f = {(a,b) : ~~} 
def discreteFunction (f ,x) :
    for i in range(0,len(f)) : # why not len(f)-1?
        if x == f[i][0] :
            return f[i][1]

    print('undefined at: '+str(x)) # debug
    return "{}"


if len(sys.argv) == 1 :
    print("defintegral [start] [end]\n\n"
          "options:\n"
          '-f   ordered pair table here\n'
          '-fun function here\n'
          '-s composite simson\n'
          '-t composite trap\n'
          '-e set error bound\n'
          '-i set iterations\n')
    exit(1)

fset = '-f'
function = '-fun'
simsonop = '-s'
trapop = '-t'
error = '-e' # calc n off of precision
iterations = '-i'

a = eval(sys.argv[1])
b = eval(sys.argv[2])
n = 5 #default
        
# the beauty of functional programming
if argexists(fset) :
    fun = eval(sys.argv[argexists(fset)+1])
    f = lambda x : discreteFunction(fun,x)

if argexists(function) :
    fun = sys.argv[argexists(function)+1]
    f = lambda x : eval(fun)

if argexists(iterations) :
    n = eval(sys.argv[argexists(iterations)+1])

elif argexists(error) :
    e = sys.argv[argexists(error)+1]
    if argexists(simsonop) :
        n = comSimErr(e)
    elif argexists(trapop) :
        n = comTrapErr(e)

if argexists(simsonop) :
    print(compositeSimson(f,a,b,n))

if argexists(trapop) :
    print(compositeTrap(f,a,b,n))

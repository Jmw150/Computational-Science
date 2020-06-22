#!/usr/bin/env python

#Various numerical analysis algorithms
# this file is a glob of all of the algorithms


# needs: a separate function to calculate error bound
#        a lambda'ing function to use this
#        an interface for easy use

# for some minor symbolic manipulation
#libraries
# libraries 
# libraries 
# libraries 
# libraries 
import sys
import sympy as sp
from numpy import *
from fractions import Fraction
import sys
import sympy as sp
from numpy import *
from fractions import Fraction
import sys
import sympy as sp
from numpy import *
from fractions import Fraction
import sys
from math import *
from fractions import Fraction
import sys
import sympy as sp
from numpy import *
from fractions import Fraction
import matplotlib.pyplot as plt #library for plotting
from numpy import *
import sys
from sympy import *
from math import * # same as using include and namespace
from sympy import *
import sys
import sys
from math import *
from fractions import Fraction
import numpy as np
from math import * # same as using include and namespace
from sympy import *
import sys
import sys
from math import *
from fractions import Fraction
import sys
from math import * 
# to save a little time fractions are not linked unless used
if argexists(fractions) :
    from fractions import Fraction
import os
import argparse # Useful for command line 
from ast import literal_eval
import sys
import glob # for regular expression matches
import tqdm # for status bar

# Spills math and sympy functions into the namespace so they 
# can be used on input
from math import * # same as using include and namespace
from sympy import *
from mpmath import fac #factorial
from sympy import *
import sys
import glob,tqdm
import numpy as np
import os
import argparse # Useful for command line 
from ast import literal_eval

from math import * # same as using include and namespace
from sympy import *
from mpmath import fac #factorial
import sys
from sympy import *
import sys
from fractions import Fraction as Frac

def factorial (x) :
    if x > 0 :
        return x * factorial(x-1)
    elif x == 0 :
        return 1
    elif x < 0 :
        return "false"

# generates taylor polynomail as a string
# formula:
#    f(x) = sum_(n=0 to inf) f^(n)(a)(x-a)^n/n!
# degree assumed to be an int, center may be a string or number
# fn == f^(n), these can be real derivatives or approximations
def genericTaylorP(degree, center, preimage) :
    terms = []
    for n in range(0,degree) :
        f = "f"+str(n)
        #+"("+str(center)+")" #not good wth indexes or (x)
        fact = str(factorial(n))
        diff = "("+str(preimage)+"-("+str(center)+"))**"+str(n)

        # nice little trick using sympy, assuming x and h used (check?)
        terms.append(str(simplify("("+f+"*"+diff+")/"+fact)))
    poly = ''
    poly += terms[0]
    for i in range(1,degree) :
        poly += " + "
        poly += terms[i]
    return str(simplify(poly)) # one last clean up

# multiple taylor polynomials around a single point
# good for making difference equations
#args:int,string,string 
#return polynomial array (sympy)
def taylorPs(smoothness, center, values) :
    polys = []
    for i in range(0,len(values)) :
        if values[i] != center :
            polys.append(poly(genericTaylorP(smoothness,center,values[i])))
#    print("polys: ")
#    for i in range(0,len(polys)):
#       print(str(polys[i])) #debug

    return polys

# pulls off the coefficients and puts them in a matrix
def matrixC(polys) :
    # make a list x list of coefficients
    coef= []
    for i in range(0,len(polys)) :
        polys[i] = poly(polys[i]) # a.coeffs() needs sympy object
        coef.append(polys[i].coeffs())

#    # remove 1 term t polynomial, it is only f(x)
#    for i in range(0,len(coef)-1) :
#        if len(coef[i]) == 1 :
#            coef.pop(i)
    return Matrix(coef)

# returns a derivative approximation formula
def differenceFormula(degree, smoothness, center, values) :

    # to allow coefficient manipulation
    x,h = symbols('x h')
    f = []
    for i in range(0,smoothness):
        f.append(Symbol("f"+str(i)))

    # both list greatest to smallest degree
    polys = taylorPs(smoothness, center, values)
    coef = matrixC(polys)

    degreepart = coef[:,degree]
    errorpart = coef[:,-1]

    # remove top degree, degree looking for, and f(x)
    # error if smoothness < 3
    indexDegree = smoothness - degree 
    coef.col_del(degree) #f to solve for
    coef.col_del(smoothness-2) # f(x)
    coef.col_del(0) #high degree

    toSolve = coef.T
    null = toSolve.nullspace()
    null = null[0]

    # this part needs the beginning parts of the equation
    # f(x+h) and the like
    # scale poly by null, add rows of poly, solve for degree
    for i in range(0,len(polys)):
        polys[i] = int(null[i])*polys[i]

    polysum = polys[0]
    for i in range(1,len(polys)):
        polysum += polys[i]

    rhsList = []
    i = 0
    j = 0
    while i < len(values):
        if values[i] == center :
            i += 1
        else :
            rhsList.append(str(-1*null[j])+"*f("+str(values[i])+")")
            i += 1
            j += 1

    rhs = rhsList[0]
    for i in range(1,len(rhsList)):
        rhs += '+'
        rhs += rhsList[i]

    # to calc bottom h value, if not symbol like x
    if type(eval(center)) == int :
        botH = 0
        for i in range(0,len(degreepart)):
            botH += Frac(values[i]) - Frac(center)
    else :
        botH = 'h'

    rhs = "("+str(-1*sum(null))+"*f("+center+")+"+rhs+")"
    rhs += "/("+str(-1*sum(null)*sum(degreepart))+'*'
    rhs += str(simplify(str(botH)+"**"+str(degree)))+")"

    return rhs

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


# function is a string, rest are numbers
# output is a sympy object
def arcLength (function, start, end) :
    x = Symbol('x')
    functionp = str(diff(function,x))
    lengthStr = str(integrate("(1 + " + functionp + "**2)**(1/2)", x))
    lengthInt = lambda x : eval(lengthStr)
    return (lengthInt(end) - lengthInt(start))


"""
Bisection method implementation

made by: Jordan Winkler
"""


# main function of the program
def bisection (f, start, stop, iterations) :
    """
    bisection [function] [start] [stop] [iterations]
 
    checks if the bisection method would work on a function of 
    a given domain
 
    impure function: 
      takes in a function and its domain
      returns True/False if the bisection method would work
      behavior changes if -f, or -d are terminal arguments
    """
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



# catch common errors before starting

    parser = argparse.ArgumentParser()
    parser.add_argument('--function',
        help='Enter smooth function to taylor expand in terms of x',
        default='sin(x)')
    parser.add_argument('--n-derivatives',
        help='number of expansions for a taylor polynomial'
        default=5)
    parser.add_argument('--center',
        help='pick point to expand on'
        default=0)
    parser.add_argument('--compute',
        help='pick number to compute with taylor polynomial'
        default=False)

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
#!/usr/bin/env python


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


ITERATION_LIMIT = 100

# initialize the matrix
A = np.array([[1,2,-2],[1,1,1],[2,2,1]])

#A = np.array([[10., -1., 2., 0.],
#              [-1., 11., -1., 3.],
#              [2., -1., 10., -1.],
#              [0., 3., -1., 8.]])
# initialize the RHS vector
b = np.array([7, 2, 5])

print("System of equations:")
for i in range(A.shape[0]):
    row = ["{0:3g}*x{1}".format(A[i, j], j + 1) for j in range(A.shape[1])]
    print("[{0}] = [{1:3g}]".format(" + ".join(row), b[i]))

x = np.zeros_like(b)
for it_count in range(1, ITERATION_LIMIT):
    x_new = np.zeros_like(x)
    print("Iteration {0}: {1}".format(it_count, x))
    for i in range(A.shape[0]):
        s1 = np.dot(A[i, :i], x_new[:i])
        s2 = np.dot(A[i, i + 1:], x[i + 1:])
        x_new[i] = (b[i] - s1 - s2) / A[i, i]
    if np.allclose(x, x_new, rtol=1e-8):
        break
    x = x_new

print("Solution: {0}".format(x))
error = np.dot(A, x) - b
print("Error: {0}".format(error))
#!/usr/bin/env python

# hermite's interpolation
# made by: Jordan Winkler
# finds the interpoling polynomial of a function 
# implementation of this function is going to just be a more
# robust interpol
# option 
#        -d   degree of interpol
#        -fn  f^(n) table here
#        -at  give back the value of P(x)

# x, f, and fp should be global
# f and fp are functions, so they should be sets of ordered pairs
# the variable "function" should be a lambda, functionStr for string
#
# last version prefered deriviatives over f data points
# this one should not overwrite any given data, (at most duplicate it)

fntable = "-f"

useDegree = "-d"
valueOnly = "-at"

xInterval = "-xin"


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


# checks if the begginging part of a string is in a list of strings
def bargexists (string, strings, start) :
    # check each word in the list
    length = len(strings)
    i = start
    while i < length :
        # check each letter

        # test if the test smaller or the user input
        depth = len(string)
        depth2 = len(strings[i])
        if depth2 < depth :
            depth = depth2

        matches = 1
        j = 0
        while j < depth :
            if strings[i][j] != string[j] :
                matches = 0
                break
            j += 1
        if matches == 1 : # if matches
            return i
        i += 1
    return 0

# long-line-itis
def getPrime (userInput, after) :
        return eval(userInput[bargexists(fntable,userInput,after)][len(fntable)])

# make function tree f[][][]
def fillf (userInput) :
    f = []
    # assume if -f then -f(number)
    k = 0 # assuming name of this function is 1
    i = 0
    while k < len(userInput) :
        if bargexists(fntable,userInput,k) :
            primes = getPrime(userInput, k)
            while i <= primes : # 0 to n
                f.append(['']) #only want to make room for max d, ah well
                i += 1
            f[primes] = eval(userInput[bargexists(fntable,userInput,k)+1])
        k += 2
    return f

userInput = sys.argv

f = fillf(userInput)
#print ("start")
#print (f)

# layer n <= 2 + layer n+1
def fillBottom (f) :
    i = len(f) - 1 - 1
    while i >= 0 :
        while len(f[i]) < len(f[i+1]) + 1 :
            f[i].append([''])
        i -= 1
    return f

f = fillBottom(f)
#print ("adjusted bottom")
#print (f)


# tree should be len(f[0]) tall
def allocPerfect (f) :
    n = len(f[0]) - 1
    i = len(f) - 1
    k = 1
    while i < n :
        f.append([])
        i += 1
    while i > 0 :
        j = len(f[i]) # j not index
        while j < k:
            f[i].append([''])
            j += 1
        i -= 1
        k += 1
    return f

f = allocPerfect(f)
#print ("made perfect tree nodes")
#print (f)

# fill down is broken as of version 2.0
# recursive perfect tree filler
# f, i, j, pair == f[i][j] = (f[i][j][0],f[i][j][1])original
def fillDown(f, i, j, pair) : 
    if i == 0 : #at leaf
        f[i][j] = pair
        return f
    else : 
        f[i][j] = pair
        # print (str(i)+ str(j) + str(f[i][j])) # debug
        f = fillDown(f,i-1,j,pair)
        f = fillDown(f,i-1,j+1,pair)
        return f

def perfectTree (f) :
    i = len(f) - 1
    while i > 0 :
        j = len(f[i]) - 1
        while j >= 0 : 
            #print("ij" + str(i) + str(j)) # debug
            # needs to be fixed
            if f[i][j][0] == f[i-1][j][0] : 
                f[i-1].insert(j+1,f[i-1][j])
                if f[i-1][-1] == [''] :
                    f[i-1].pop() 
                if f[i][j] != f[i][-1] :
                    f[i].insert(j+1,[''])
                #f = fillDown(f,i,j,(f[i][j][0],f[i][j][1]))
            j -= 1
        i -= 1
    return f

f = perfectTree(f)
f = allocPerfect(f)
# debug
#print ("filled out tree") 
#print (f)

# divided difference notation
def F(i,j) :
   #print ("ij proper:" + str(i) + str(j))
   #print ("ij:" + str(j) + str(i-j))
   if f[j][i-j] != [''] :
       return f[j][i-j][1] #index F(i,j) is reversed...
   else :
       return (F(i,j-1) - F(i-1,j-1))/float(f[0][i][0] - f[0][i-j][0])

def makeNewtonPoly (degree=len(f[0])) :
    if degree < len(f[0])-1 : # valid if want less, not more
        n = degree
    else :
        n = len(f[0]) - 1
    P = str(F(0,0))

    i = 1
    while i <= n :
        #print ("i:" + str(i)) #debug
        P = P + "+(" + str(F(i,i)) + ")"
        i += 1
        
        j = 0
        while j < (i-1) : #did one more time than asked
            P = P + "*(x-(" + str(f[0][j][0]) + "))"
            j += 1

    return P

if argexists(useDegree) :
    degree = eval(sys.argv[argexists(useDegree)+1])
else :
    degree = len(f[0])

P = makeNewtonPoly(degree)
newPoly = lambda x : eval(P)


if argexists(valueOnly) :
    xval = eval(sys.argv[argexists(valueOnly)+1])
    print(newPoly(xval))
else :
    print (P)

# not my code, just used an example
def horner(x0, *a): #does this even work?
    '''
        Horner's method is an algorithm to calculate a polynomial at
        f(x0) and f'(x0)

        x0 - The value to avaluate
        a - An array of the coefficients

        The degree is the polynomial is set equal to the number of coefficients
    '''
    n = len(a)

    y = a[0]
    z = a[0]
    for j in range(1, n - 1):
        y = x0 * y + a[j]
        z = x0 * z + y

    y = x0 * y + a[-1]

    print('P(x0) =', y)
    print('P\'(x0) =', z)

# not my code end


def hornersMethod(x_0, *a) : #broken?
    n = len(a)
    
    y = a[0]
    z = a[0]
    
    i = 0
    while (i < n) :
        y = x_0 * y + a[i]
        z = x_0 * z + y
        i += 1

    y = x_0 * y + a[-1]

    print("P(x_0)  = " + str(y))
    print("P\'(x_0) = " + str(z))

def poly_horner(A, x): #just an example, not my code
    p = A[-1]
    i = len(A) - 2
    while i >= 0:
        p = p * x + A[i]
        i -= 1
    return p


poly = ( 1, 0, 0) #x^2

x = 2

print (poly_horner(poly, x))
#!/usr/bin/env python


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
#!/usr/bin/env python


def norm1 (v) :
    norm = 0
    for i in range(0,len(v)):
        norm += abs(v[i])
    return norm

def norm (v) :
    norm = 0
    for i in range(0,len(v)):
        norm += v[i]**2
    norm**0.5
    return norm

def norminf (v) :
    for i in range(0,len(v)):
        v[i] = abs(v[i])
    return max(v)

def norminfM (A) :
    a = A[0][:] # get length
    for i in range(0,len(a)): #clean it
        a[i] = 0
    for i in range(0,len(a)):
        for j in range(0,len(a)):
            a[i] += abs(A[i][j])
    return max(a)

def norm1M (A) :
    a = A[0][:] # get length
    for i in range(0,len(a)): #clean it
        a[i] = 0
    for i in range(0,len(a)):
        for j in range(0,len(a)):
            a[i] += abs(A[j][i])
    return max(a)

# jacobi iterative
# input: n = len(A[0]), A, b, x=guess, tolerance or iterations
def jacobi (A, x, b, tolerance, iterations) :
    xo = x
    k = 1 # step 1
    n = len(A[0]) 
    while (k <= iterations) : # step 2 
        for i in range(0,n) : # step 3 
            p1 = 0
            for j in range(0,n) :
                if (j != i) : 
                    p1 += A[i][j]*xo[j]
            x[i] = (1/float(A[i][i]) * (-p1 + b[i]))
      
        xdiff = [0] * n # step 4
        for i in range(0,n) :
            xdiff[i] = x[i] - xo[i]
        if (norm(xdiff) < tolerance and k > 1) :
            return (k, x)

        k += 1 # step 5
        xo = x # step 6
    return "err: " + str(x) # step 7

# n = len(A[0]), A, b, xo, tol, iterations
def gaussSeidel (A, x, b, tolerance, iterations) : # wrong
    n = len(A[0])
    xo = x
    k = 1 # step 1
    while (k <= iterations) : # step 2
        for i in range(0,n) : # step 3
            p1 = 0
            p2 = 0
            for j in range(0,i) : # range does not reach end
                p1 += A[i][j]*x[j]
            for j in range(i+1, n) :
                p2 += A[i][j]*xo[j]
            x[i] = 1/float(A[i][i])*(-p1-p2+b[i])

        xd = [0] * n # step 4
        for i in range(0,n) :
            xd[i] = x[i] - xo[i]
        if (norm(xd) < tolerance) :
            return (k, x)

        k += 1 # step 5
        xo = x # step 6
    return "err: " + str(x) # step 7
               
A = [[1.0,2.0,-2.0],[1.0,1.0,1.0],[2.0,2.0,1.0]]
b = [7.0,2.0,5.0]
x = [0,0,0]
iteration = 1
tolerance = 10**(-5)

A2 = [[1,1,2],[2,1,2],[-2,1,1]]
A1 = [[3,-1,1],[2,5,2],[-1,-1,3]]
b1 = [0,10,-6]
A2 = [[1,2,-2],[1,1,1],[2,2,1]]
b2 = [7,2,5]
A3 = [[2,-1,1],[2,2,2],[-1,-1,2]]
b3 = [-1,4,-5]

# actual answer [3, 1/2, -3/2]
# 1.
#print ("jac: " +str(jacobi(A, x, b, tolerance, 1)))
#print ("jac: " +str(jacobi(A, x, b, tolerance, 2)))
# 2.
#print ("gSeid: " +str(gaussSeidel(A, x, b, tolerance, 1)))
#print ("gSeid: " +str(gaussSeidel(A, x, b, tolerance, 2)))

print (jacobi(A1,x,b1,tolerance,50))
print (jacobi(A2,x,b2,tolerance,50))
print (jacobi(A3,x,b3,tolerance,50))
#
print (gaussSeidel(A1,x,b1,tolerance,50))
print (gaussSeidel(A2,x,b2,tolerance,50))
print (gaussSeidel(A3,x,b3,tolerance,50))


#!/usr/bin/env python

# lagrange interpolation
# made by: Jordan Winkler
# makes a polynomial that interpolates data points in 2d
# This program is weak against injection attacks.

# x:list, y:list, degree:integer
def lagrangeInterpol (x, y, degree) :
    # L_i = product_{k != i, k=0 to n} (x-x_k)/(x_i - x_k)
    L = []
    i = 0
    k = 0
    while (i < degree) :

        i += 1
    # P(x) = sum_{i=0 to n} y_i L_i
    

def getL (x, degree)
    
#!/usr/bin/env python


if len(sys.argv) != 4 :
    print ("plot [function] [start-xy] [end-xy]")

start = eval(sys.argv[2])
end = eval(sys.argv[3])

#arg start domain, end domain, precision of computation
x = arange(start, end, 0.01) 

f = sys.argv[1]

y = eval(f) 

#pass t object and s function to plot. Stored at run time?
#auto adjusts to include vertical values
plt.plot(x, y) 

# set max and min values
axes = plt.gca()
axes.set_xlim([start,end])
axes.set_ylim([start-0.01,end+0.01])

plt.xlabel('x')
plt.ylabel('y')
plt.title(f)
plt.grid(True)
#plt.savefig("graph.png")
plt.show() #call a program to display this grid

#!/usr/bin/env python

# Muller's method
# made by: Jordan Winkler
# finds an approximation of a root using Muller's method
# This program is weak against injection attacks.
# option -f for fractions, default floating point
#        -d for debugging or verbose mode
#        -i [number] for iterations of function
#        -p [number] for precision

# could have been defined recursively or with a do while loop
# f is a function, the rest are numbers
# returns number on success, string on failure
def mullersMethod (f, p0, p1, p2, iterations, tolerance) :
    # set up parabola interpolation
    h1 = p1 - p0
    h2 = p2 - p1
    s1 = (f(p1) - f(p0))/float(h1)
    s2 = (f(p2) - f(p1))/float(h2)
    d  = (s2 - s1)/float((h2 + h1))
    i = 3

    while i <= iterations :
        # calc base
        b = s2 + h2*d
        D = ((b**2 - 4*f(p2)*d)+0j)**0.5 #convert to complex 
        
        # pick the larger base
        if abs(b-D) < abs(b+D) :
            E = b + D
        else :
            E = b - D

        # calc rest of equation
        h = -2*f(p2)/E
        p = p2 + h
        
        # if within tolerance
        if abs(h) < tolerance :
            return p

        # update values to go again
        p0 = p1
        p1 = p2
        p2 = p
        h1 = p1 - p0
        h2 = p2 - p1
        s1 = (f(p1) - f(p1))/(h1)
        s2 = (f(p2) - f(p1))/(h2)
        d  = (s2 - s1)/(h2 + h1)
        i += 1
    
    return " error "
        



fstring = sys.argv[1]
f = lambda x : eval(fstring)
p0 = eval(sys.argv[2])
p1 = eval(sys.argv[3])
p2 = eval(sys.argv[4])
iterat = 50
tol = 0.0001

print (mullersMethod(f, p0, p1, p2, iterat, tol))
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




if __name__ == '__main__':

    # This code needs at least the function to run
    
    if len(sys.argv) == 1 :
        print ("taylorPolynomial [function] (n-derivatives) (centered on)")
        exit()

    parser = argparse.ArgumentParser()
    parser.add_argument('function', metavar='f', type=str,
        help='Enter smooth function to taylor expand in terms of x')
    parser.add_argument('--function',
        help='Enter smooth function to taylor expand in terms of x',
        default=sys.argv[1]) # assume it is in the first spot if not specified
    parser.add_argument('--n',
        help='number of expansions for a taylor polynomial',
        default='5')
    parser.add_argument('--center',
        help='pick point to expand on',
        default='0')
    parser.add_argument('--eval',
        help='pick number to compute with taylor polynomial',
        default='False')

    args = parser.parse_args()
    

    def Taylor_polynomial(f,a,n):
        """
        A taylor polynomial off of the series equation
        sum_{n=0}^\inf f^(n)(a)/n! * (x-a)^n

        takes: f,a,n
        f is some differentiable function
        a is some real number
        n is some positive integer
        """
        x = Symbol('x')

        # make some derivatives first, since f^(n) is recursive
        derivative = []
        derivative.append(f)
        index = 0
        while index <= n :
            derivative.append(str(diff(derivative[index],x)))
            index += 1
        #print(derivative) #debug

        taylor = ''
        index = 0
        while index <= n :
            taylor += str((lambda x : eval(derivative[index]))(a)) + '/'+str(fac(index))+' * '+str((x-a)**index)+' + '
            index += 1
        taylor = str(eval(taylor[:-3]))

        return taylor


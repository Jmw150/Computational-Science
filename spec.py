#!/usr/bin/env python

#Various numerical analysis algorithms
# this file is a glob of all of the algorithms

import sys
import os
import sympy as sp
import numpy as np
import math
from mpmath import fac #factorial
from fractions import Fraction
from fractions import Fraction as Frac
import matplotlib.pyplot as plt #library for plotting
import argparse # Useful for command line 
from ast import literal_eval
import glob # for regular expression matches
import tqdm # for status bar


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

# function is a string, rest are numbers
# output is a sympy object
def arcLength (function, start, end) :
    x = Symbol('x')
    functionp = str(diff(function,x))
    lengthStr = str(integrate("(1 + " + functionp + "**2)**(1/2)", x))
    lengthInt = lambda x : eval(lengthStr)
    return (lengthInt(end) - lengthInt(start))


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

# function as a set of ordered pairs
# set f = {(a,b) : ~~} 
def discreteFunction (f ,x) :
    for i in range(0,len(f)) : # why not len(f)-1?
        if x == f[i][0] :
            return f[i][1]

    print('undefined at: '+str(x)) # debug
    return "{}"


# hermite's interpolation
# made by: Jordan Winkler
# finds the interpoling polynomial of a function 
# implementation of this function is going to just be a more
# robust interpol
# option 
#        -d   degree of interpol
#        -fn  f^(n) table here
#        -at  give back the value of P(x)


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

# layer n <= 2 + layer n+1
def fillBottom (f) :
    i = len(f) - 1 - 1
    while i >= 0 :
        while len(f[i]) < len(f[i+1]) + 1 :
            f[i].append([''])
        i -= 1
    return f

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

# jacobiXiterative
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
               
# lagrange interpolation
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


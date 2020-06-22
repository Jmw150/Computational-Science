#!/usr/bin/env python

# needs: a separate function to calculate error bound
#        a lambda'ing function to use this
#        an interface for easy use

# for some minor symbolic manipulation
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

# terminal arguments
derive = "-d"
cent = "-c"
val = "-v"
smooth = "-s"

# data
if len(sys.argv) != 4 :
    print("approxderivative [derivative] [center] [values] [smoothness]")
    exit(1)
else : 
    derivative = eval(sys.argv[1])
    center = sys.argv[2]
    values = eval(sys.argv[3]) # so bash does not freak out
    smoothness = eval(sys.argv[4])

    a = differenceFormula(derivative, smoothness, center, values)
    print(a)


#derivative = 1 # derivative to get back
#center = "x" 
#values = ["x-h","x","x+h","x+2*h"] 
#smoothness = 5

#minH = genericTaylorP(4,"x","x-h")
#plusH = genericTaylorP(4,"x","x+h")
#print (plusH)
#print (minH)
#print (simplify('(('+plusH+')-('+minH+'))/(2*h)'))
#
#print (genericTaylorP(3,"x","x"))

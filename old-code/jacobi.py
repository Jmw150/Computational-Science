#!/usr/bin/env python

from sympy import *

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



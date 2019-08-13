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

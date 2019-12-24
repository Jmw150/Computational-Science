#!/usr/bin/env python

#

import glob,tqdm
import numpy as np
import os
import argparse # Useful for command line 
from ast import literal_eval

from math import * # same as using include and namespace
from sympy import *
from mpmath import fac #factorial
import sys


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

    if args.eval != False :
        print(
            (lambda x : eval(
                Taylor_polynomial(args.function, 
                                  eval(args.center), 
                                  eval(args.n))))(eval(args.eval)))
    else :
        print(Taylor_polynomial(args.function, 
                                eval(args.center), 
                                eval(args.n)))

    



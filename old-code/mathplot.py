#!/usr/bin/env python

import matplotlib.pyplot as plt #library for plotting
from numpy import *
import sys

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


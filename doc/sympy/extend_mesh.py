#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
from sympy.matrices import *
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False, no_global=True)
from sympy.printing import cxxcode

s,x,alpha = symbols('s,x,alpha',positive=True)

y = d*x**3 + c*x**2 + b*x + a
ydiff = diff(y,x)
print("ansatz: y = {}\n        y' = {}".format(y, ydiff))

params = solve([y.subs(x,1-s)-(1-s), y.subs(x,1-s/2)-(1-s+alpha*s/2), y.subs(x,1)-1, ydiff.subs(x,1-s)-1], [a,b,c,d])
print("\nparams: {}".format(params))

y = y.subs(params)
ydiff = ydiff.subs(params)
print("\nresult: y = {}".format(y))
print("\n{}".format(cxxcode(y, standard='C++11')))

from matplotlib import pyplot as plt
import numpy as np
svalue = 0.95
avalue = 0.7
xlist = np.linspace(1-svalue, 1, 20)
ylist = [y.subs([(alpha,avalue),(s,svalue),(x,xvalue)]).evalf() for xvalue in xlist]
ydifflist = [ydiff.subs([(alpha,avalue),(s,svalue),(x,xvalue)]).evalf() for xvalue in xlist]

print("numerical function: {}".format(y.subs([(s,svalue),(alpha,avalue)])))

# plot curves
# set global parameters for font sizes
plt.rcParams.update({'font.size': 14})
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 8
# colors: #800000,   #ffcc00,     #ff7f0e,   #d62728
fig,axs = plt.subplots(2,1,sharex=True)
axs[1].set_xlabel('r')
axs[0].grid(which='both')
axs[0].plot([0,1-svalue], [1,1], "-", color="#800000")
axs[0].plot(xlist, ydifflist, "-", color="#800000", label="$y'$")
axs[0].legend()

axs[1].plot(xlist, ylist, "-", color="#ff7f0e", label="$y$")
axs[1].plot([1-svalue/2,1-svalue/2],[1-svalue,1-svalue/2], "-", color="#ffcc00")
axs[1].plot([0,1-svalue],[0,1-svalue], "-", color="#ff7f0e")
axs[1].plot([1-svalue,1],[1-svalue,1], ":", color="#ff7f0e")
axs[1].grid(which='both')
axs[1].plot([1-svalue/2], [1-svalue+avalue*svalue/2], "o", color="#ffcc00")
axs[1].plot([1-svalue], [1-svalue], "o", color="#ff7f0e")
axs[1].plot([1], [1], "o", color="#ff7f0e")
axs[1].plot([1-svalue/2,1], [1-svalue/2,1], "o", color="#ff7f0e")
axs[1].plot([1-svalue/2], [1-svalue], "o", color="#ff7f0e")
axs[1].legend()

plt.savefig("extend_mesh_plot.pdf")
plt.show()

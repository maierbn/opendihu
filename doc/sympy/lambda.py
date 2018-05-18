#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
from sympy.matrices import *
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False, no_global=True)


c0, lambdaValue,tmax,that,lz,ly = symbols('c0,lambdaValue,tmax,that,lz,ly',positive=True)

if False:
  print "(3D) total load constant (traction in reference configuration), direct solution"
  print solve(2*c0*lambdaValue**3 - tmax/(lz*lz)*lambdaValue**2 - 2*c0, lambdaValue, simplify=True, positive=True)[2]

  print "(3D) total load constant (traction in reference configuration), penalty formulation"
  print solve(2*c0*lz*lz*(lambdaValue**3 - 1) - tmax*lambdaValue**2, lambdaValue, simplify=True, positive=True)[2]

  print "(3D) traction constant"
  print solve(2.*c0*lambdaValue**3 - that*lambdaValue - 2.*c0, lambdaValue, simplify=True, positive=True)[0]

print "(2D) total load constant (traction in reference configuration), penalty formulation"
term_penalty = solve(2.*c0*ly*(lambdaValue**4. - 1.) - tmax*lambdaValue**3., lambdaValue, simplify=True, positive=True)[3]

print term_penalty
f = term_penalty.subs([(ly,1.0), (c0,1.0)])

from matplotlib import pyplot as plt
import numpy as np
xlist = np.linspace(0.1, 2.0, 20)
ylist = [f.subs(tmax,x).evalf() for x in xlist]

print xlist
print ylist

plt.xlabel('traction')
plt.ylabel('lambda')
plt.plot(xlist, ylist)
plt.show()

print term_penalty.subs([(ly,1.0), (c0,1.0), (tmax,0.5)]).evalf()

print "(2D) total load constant (traction in reference configuration),  direct solution"
term_direct = solve(2*c0*lambdaValue**4 - tmax/ly*lambdaValue**3 - 2*c0, lambdaValue, simplify=True, positive=True)[3]


print term_direct
print term_direct.subs([(ly,1.0), (c0,1.0), (tmax,0.5)]).evalf()

print "alternative parametrization"

lx_settings = 1.5
ly_settings = 0.6
tmax_settings = 2.2  
analytic_lambda = term_direct.subs([(ly,ly_settings), (c0,1.0), (tmax,tmax_settings)]).evalf()

print "analytic_lambda:",analytic_lambda

# displacements solution for 2x2 quadratic elements
reference_solution = [
  # bottom row
  np.array([0.00,                             0.0]), 
  np.array([0.25*analytic_lambda*lx_settings, 0.0]),
  np.array([0.50*analytic_lambda*lx_settings, 0.0]),
  np.array([0.75*analytic_lambda*lx_settings, 0.0]),
  np.array([1.00*analytic_lambda*lx_settings, 0.0]),
  
  np.array([0.00,                             0.25/analytic_lambda*ly_settings]), 
  np.array([0.25*analytic_lambda*lx_settings, 0.25/analytic_lambda*ly_settings]),
  np.array([0.50*analytic_lambda*lx_settings, 0.25/analytic_lambda*ly_settings]),
  np.array([0.70*analytic_lambda*lx_settings, 0.25/analytic_lambda*ly_settings]),
  np.array([1.00*analytic_lambda*lx_settings, 0.25/analytic_lambda*ly_settings]),
  
  np.array([0.00,                             0.5/analytic_lambda*ly_settings]), 
  np.array([0.25*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings]),
  np.array([0.50*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings]),
  np.array([0.70*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings]),
  np.array([1.00*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings]),
  
  np.array([0.00,                             0.75/analytic_lambda*ly_settings]),
  np.array([0.25*analytic_lambda*lx_settings, 0.75/analytic_lambda*ly_settings]),
  np.array([0.50*analytic_lambda*lx_settings, 0.75/analytic_lambda*ly_settings]),
  np.array([0.70*analytic_lambda*lx_settings, 0.75/analytic_lambda*ly_settings]),
  np.array([1.00*analytic_lambda*lx_settings, 0.75/analytic_lambda*ly_settings]),
  # top row
  np.array([0.00,                             1.00/analytic_lambda*ly_settings]),
  np.array([0.25*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings]),
  np.array([0.50*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings]),
  np.array([0.70*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings]),
  np.array([1.00*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings]),
]

# displacements solution for 1x2 quadratic elements
reference_solution = [
  # bottom row
  np.array([0.00,                             0.0]), 
  np.array([0.25*analytic_lambda*lx_settings, 0.0]),
  np.array([0.50*analytic_lambda*lx_settings, 0.0]),
  np.array([0.75*analytic_lambda*lx_settings, 0.0]),
  np.array([1.00*analytic_lambda*lx_settings, 0.0]),
  # center row
  np.array([0.00,                             0.5/analytic_lambda*ly_settings]), 
  np.array([0.25*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings]),
  np.array([0.50*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings]),
  np.array([0.70*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings]),
  np.array([1.00*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings]),
  
  # top row
  np.array([0.00,                             1.00/analytic_lambda*ly_settings]),
  np.array([0.25*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings]),
  np.array([0.50*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings]),
  np.array([0.70*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings]),
  np.array([1.00*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings]),
]
print reference_solution












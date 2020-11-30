#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
from sympy.matrices import *
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False, no_global=True)
from sympy.printing import cxxcode
import sys

#  ^   6  
#  |   3  8
# xi2  0__1__2
#  |   --> xi1

# compute ansatz functions
xi1,xi2 = symbols('xi1,xi2', real=True)
c00,c10,c01,c11,c20,c02 = symbols('c00,c10,c01,c11,c20,c02', real=True)

ansatz = c00 + c10*xi1 + c01*xi2 + c11*xi1*xi2 + c20*xi1*xi1 + c02*xi2*xi2
print("ansatz: {}".format(ansatz))

nodal_coordinates = [(0,0), (1/2,0), (1,0), (0,1/2), (0,1), (1/2,1/2)]
ansatz_functions = []

for i in range(6):
  conditions = []
  
  # loop over nodes
  for k,(xi1_coordinate,xi2_coordinate) in enumerate(nodal_coordinates):
    value = 0
    if i==k:
      value = 1
    
    conditions.append(ansatz.subs({xi1: xi1_coordinate, xi2: xi2_coordinate}) - value)

  constants = solve(conditions, [c00,c10,c01,c11,c20,c02])
  constants = simplify(constants)
  #print("constants: {}".format(constants))
  
  ansatz_function = ansatz.subs(constants)
  ansatz_functions.append(ansatz_function)
  
  print("ansatz function {}: {}".format(i,factor(ansatz_function)))

# compute jacobian
xp1,x01,x11,x21,x31,x41,x51 = symbols('xp1,x01,x11,x21,x31,x41,x51', real=True)
xp2,x02,x12,x22,x32,x42,x52 = symbols('xp2,x02,x12,x22,x32,x42,x52', real=True)

Phi1 = ansatz_functions[0]*x01 + ansatz_functions[1]*x11 + ansatz_functions[2]*x21 + ansatz_functions[3]*x31 + ansatz_functions[4]*x41 + ansatz_functions[5]*x51
Phi2 = ansatz_functions[0]*x02 + ansatz_functions[1]*x12 + ansatz_functions[2]*x22 + ansatz_functions[3]*x32 + ansatz_functions[4]*x42 + ansatz_functions[5]*x52

jacobian = [[diff(Phi1,xi1), diff(Phi1,xi2)],
            [diff(Phi2,xi1), diff(Phi2,xi2)]]

print("jacobian J_11: {}".format(simplify(jacobian[0][0])))
print("jacobian J_12: {}".format(simplify(jacobian[1][0])))
print("jacobian J_21: {}".format(simplify(jacobian[0][1])))
print("jacobian J_22: {}".format(simplify(jacobian[1][1])))
print("")
print("factorized:")
print("jacobian J_11: {}".format(factor(jacobian[0][0])))
print("jacobian J_12: {}".format(factor(jacobian[1][0])))
print("jacobian J_21: {}".format(factor(jacobian[0][1])))
print("jacobian J_22: {}".format(factor(jacobian[1][1])))

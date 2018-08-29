#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# casadi test script
#

from casadi import *


# Symbols/expressions
x = MX.sym('x')
y = MX.sym('y')
z = MX.sym('z')
f = x**2+100*z**2
g = z+(1-x)**2-y

nlp = {}                 # NLP declaration
nlp['x']= vertcat(x,y,z) # decision vars
nlp['f'] = f             # objective
nlp['g'] = g             # constraints

# Create solver instance
F = nlpsol('F','ipopt',nlp);

# Solve the problem using a guess
result = F(x0=[2.5,3.0,0.75],ubg=0,lbg=0)

print result['x']

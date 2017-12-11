#!/usr/bin/pxi2thon
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
  

xi1,xi2 = symbols('xi1 xi2')
Phi1, Phi2, Phi3 = symbols('Phi1 Phi2 Phi3')
# 3 4
# 1 2
x11, x12, x13 = symbols('x11 x12, x13')
x21, x22, x23 = symbols('x21 x22, x23')
x31, x32, x33 = symbols('x31 x32, x33')
x41, x42, x43 = symbols('x41 x42, x43')

# arbitrary 2D case

# define Phi(xi)
Phi1 = (1-xi1)*(1-xi2)*x11 + xi1*(1-xi2)*x21 + (1-xi1)*xi2*x31 + xi1*xi2*x41
Phi2 = (1-xi1)*(1-xi2)*x12 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x32 + xi1*xi2*x42
Phi3 = (1-xi1)*(1-xi2)*x13 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x33 + xi1*xi2*x43

# compute J2
dPhi1_dxi1 = diff(Phi1, xi1)
dPhi2_dxi1 = diff(Phi2, xi1)
dPhi3_dxi1 = diff(Phi3, xi1)
dPhi1_dxi2 = diff(Phi1, xi2)
dPhi2_dxi2 = diff(Phi2, xi2)
dPhi3_dxi2 = diff(Phi3, xi2)

dPhi_dxi1_squared = dPhi1_dxi1**2 + dPhi2_dxi1**2 + dPhi3_dxi1**2
dPhi_dxi2_squared = dPhi1_dxi2**2 + dPhi2_dxi2**2 + dPhi3_dxi2**2
dPhi_dxi1_dxi2 = dPhi1_dxi1*dPhi1_dxi2 + dPhi2_dxi1*dPhi2_dxi2 + dPhi3_dxi1*dPhi3_dxi2

gram_det = dPhi_dxi1_squared*dPhi_dxi2_squared - dPhi_dxi1_dxi2**2
j = sqrt(gram_det)

print "arbitrary J_2:", j

# rectangular cartesian 2D
l1, l2 = symbols('l1 l2')
Phi1 = x11 + xi1*l1
Phi2 = x12 + xi2*l2
Phi3 = 0

# compute J2
dPhi1_dxi1 = diff(Phi1, xi1)
dPhi2_dxi1 = diff(Phi2, xi1)
dPhi3_dxi1 = diff(Phi3, xi1)
dPhi1_dxi2 = diff(Phi1, xi2)
dPhi2_dxi2 = diff(Phi2, xi2)
dPhi3_dxi2 = diff(Phi3, xi2)

dPhi_dxi1_squared = dPhi1_dxi1**2 + dPhi2_dxi1**2 + dPhi3_dxi1**2
dPhi_dxi2_squared = dPhi1_dxi2**2 + dPhi2_dxi2**2 + dPhi3_dxi2**2
dPhi_dxi1_dxi2 = dPhi1_dxi1*dPhi1_dxi2 + dPhi2_dxi1*dPhi2_dxi2 + dPhi3_dxi1*dPhi3_dxi2

gram_det = dPhi_dxi1_squared*dPhi_dxi2_squared - dPhi_dxi1_dxi2**2
j = sqrt(gram_det)

print "RC J_2:", simplify(j)


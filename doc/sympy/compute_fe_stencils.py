#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *

# compute stencil for FE Δ

xi,xi1,xi2,xi3 = symbols('xi xi1 xi2 xi3')

##########################################################################
### 1D element contributions for linear Lagrange ansatz function (phi) ###
def phi0(xi):
  return (1-xi)
  
def phi1(xi):
  return xi
  
print "1D stencil"
print "=========="

# left (center) node
integrand = diff(phi0(xi),xi)**2
entry0 = integrate(-integrand,(xi,0,1))

# right node
integrand = diff(phi0(xi),xi)*diff(phi1(xi),xi)
entry1 = integrate(-integrand,(xi,0,1))

print "element contribution:"
print "[",entry0,entry1,"]"

print "nodal stencil:"
print "[",entry1,2*entry0,entry1,"]"

########################################################
### 1D element contributions for cubic Hermite (phi) ###

def phi0(xi):
  return 1 - 3*xi*xi + 2*xi*xi*xi

def phi1(xi):
  return xi * (xi-1) * (xi-1)      # xi*xi*xi - 2*xi*xi + xi

def phi2(xi):
  return xi*xi * (3 - 2*xi)        # -2*xi*xi*xi + 3*xi*xi

def phi3(xi):
  return xi*xi * (xi-1)              # xi*xi*xi - xi*xi

print "1D stencil Hermite"
print "=========="

integrand = diff(phi0(xi),xi) * diff(phi2(xi),xi)
entry02 = integrate(-integrand,(xi,0,1))

integrand = diff(phi1(xi),xi) * diff(phi2(xi),xi)
entry12 = integrate(-integrand,(xi,0,1))

integrand = diff(phi2(xi),xi) * diff(phi2(xi),xi)
entry22 = integrate(-integrand,(xi,0,1))

integrand = diff(phi3(xi),xi) * diff(phi2(xi),xi)
entry32 = integrate(-integrand,(xi,0,1))


integrand = diff(phi0(xi),xi) * diff(phi3(xi),xi)
entry03 = integrate(-integrand,(xi,0,1))

integrand = diff(phi1(xi),xi) * diff(phi3(xi),xi)
entry13 = integrate(-integrand,(xi,0,1))

integrand = diff(phi2(xi),xi) * diff(phi3(xi),xi)
entry23 = integrate(-integrand,(xi,0,1))

integrand = diff(phi3(xi),xi) * diff(phi3(xi),xi)
entry33 = integrate(-integrand,(xi,0,1))

print "entry 02: ", entry02, "=", float(entry02)
print "entry 12: ", entry12, "=", float(entry12)
print "entry 22: ", entry22, "=", float(entry22)
print "entry 32: ", entry32, "=", float(entry32)
print ""
print "entry 03: ", entry03, "=", float(entry03)
print "entry 13: ", entry13, "=", float(entry13)
print "entry 23: ", entry23, "=", float(entry23)
print "entry 33: ", entry33, "=", float(entry33)

# linear line in geometry, given by Hermite ansatz functions
c0,c1,c2,c3 = symbols('c0,c1,c2,c3')
xa,xb = symbols('xa,xb')

def x1(t):
  return (xb-xa)*t + xa
def x2(t):
  return c0*phi0(t) + c1*phi1(t) + c2*phi2(t) + c3*phi3(t)

print "x(t) =",factor(x2(t),t)
p1 = Poly(x1(t), t)
coeffs1 = p1.all_coeffs()
p2 = Poly(x2(t), t)
coeffs2 = p2.all_coeffs()
print "coefficients of x1: ", coeffs1
print "coefficients of x2: ", coeffs2

print "solve for c0,c1,c2,c3"
result = solve([coeffs2[0], coeffs2[1], coeffs2[2] - coeffs1[0], coeffs2[3] - coeffs1[1]],[c0,c1,c2,c3])
print result
print "test:"
p2 = x2(t)
print p2
print ""
print factor(p2.subs(result),t)
print ""
print "conclusion: to represent a 1D linear geometry in an element between xa and xb, use the following nodal dofs:"
print "node0: dof0 (c0):",result[c0]
print "node0: dof1 (c1):",result[c1]
print "node1: dof0 (c2):",result[c2]
print "node1: dof1 (c3):",result[c3]
  
##########################################################################
### 2D element contributions for linear Lagrange ansatz function (phi) ###

# 2 3
# 0 1

def phi0(xi1,xi2):
  return (1-xi1)*(1-xi2)
  
def phi1(xi1,xi2):
  return xi1*(1-xi2)

def phi2(xi1,xi2):
  return (1-xi1)*xi2

def phi3(xi1,xi2):
  return xi1*xi2
  
print "2D stencil"
print "=========="

# lower left node
integrand = diff(phi0(xi1,xi2),xi1)**2 + diff(phi0(xi1,xi2),xi2)**2
entry00 = integrate(-integrand,(xi1,0,1),(xi2,0,1))      # -2/3

# lower right node
integrand = diff(phi0(xi1,xi2),xi1)*diff(phi1(xi1,xi2),xi1) + diff(phi0(xi1,xi2),xi2)*diff(phi1(xi1,xi2),xi2)
entry01 = integrate(-integrand,(xi1,0,1),(xi2,0,1))      # 1/6

# upper left node
integrand = diff(phi0(xi1,xi2),xi1)*diff(phi2(xi1,xi2),xi1) + diff(phi0(xi1,xi2),xi2)*diff(phi2(xi1,xi2),xi2)
entry10 = integrate(-integrand,(xi1,0,1),(xi2,0,1))      # 1/6

# upper right node
integrand = diff(phi0(xi1,xi2),xi1)*diff(phi3(xi1,xi2),xi1) + diff(phi0(xi1,xi2),xi2)*diff(phi3(xi1,xi2),xi2)
entry11 = integrate(-integrand,(xi1,0,1),(xi2,0,1))      # 1/3

print "element contribution:"
print "[",entry10,entry11,"]"
print "[",entry00,entry01,"]"

print "nodal stencil:"
print "[",entry11,2*entry10,entry11,"]"
print "[",2*entry01,4*entry00,2*entry01,"]"
print "[",entry11,2*entry10,entry11,"]"

# the element contribution stencil is
# [   1/6   1/3 ]
# [ _-2/3_  1/6 ]

# the full stencil is therefore
# [ 1/3   1/3   1/3 ]       [ 1   1   1 ]
# [ 1/3 _-8/3_  1/3 ] = 1/3 [ 1 _-8_  1 ]
# [ 1/3   1/3   1/3 ]       [ 1   1   1 ]

## analytically compute |det(DPhi(xi))|
x11,x12,x21,x22,x31,x32,x41,x42 = symbols('x11,x12,x21,x22,x31,x32,x41,x42')

def Phi(xi1,xi2):
  p1 = (1-xi1)*(1-xi2)*x11 + xi1*(1-xi2)*x21 + (1-xi1)*xi2*x31 + xi1*xi2*x41
  p2 = (1-xi1)*(1-xi2)*x12 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x32 + xi1*xi2*x42
  
  return [p1,p2]

matrix = Matrix([[diff(Phi(xi1,xi2)[0],xi1), diff(Phi(xi1,xi2)[0],xi2)],[diff(Phi(xi1,xi2)[1],xi1), diff(Phi(xi1,xi2)[1],xi2)]])
term = det(matrix)
print "det(DPhi(xi)):"
print collect(term,(xi1,xi2))

# lower left node
integrand = (diff(phi0(xi1,xi2),xi1)**2 + diff(phi0(xi1,xi2),xi2)**2)*term
print "00", factor(integrate(-integrand,(xi1,0,1),(xi2,0,1)))

##############################################################
# linear area in geometry, given by Hermite ansatz functions #

def phi0(xi):
  return 1 - 3*xi*xi + 2*xi*xi*xi

def phi1(xi):
  return xi * (xi-1) * (xi-1)      # xi*xi*xi - 2*xi*xi + xi

def phi2(xi):
  return xi*xi * (3 - 2*xi)        # -2*xi*xi*xi + 3*xi*xi

def phi3(xi):
  return xi*xi * (xi-1)              # xi*xi*xi - xi*xi

c00,c01,c02,c03,c10,c11,c12,c13,c20,c21,c22,c23,c30,c31,c32,c33 = symbols('c00,c01,c02,c03,c10,c11,c12,c13,c20,c21,c22,c23,c30,c31,c32,c33')
x0,x1,x2,x3 = symbols('x0,x1,x2,x3')

def x_1(xi1,xi2):
  return (1-xi1)*(1-xi2)*x0 + xi1*(1-xi2)*x1 + (1-xi1)*xi2*x2 + xi1*xi2*x3
  
def x_2(xi1,xi2):
  return (c00*phi0(xi1)*phi0(xi2) + c01*phi0(xi1)*phi1(xi2) + c02*phi1(xi1)*phi0(xi2) + c03*phi1(xi1)*phi1(xi2))*x0 \
    + (c10*phi2(xi1)*phi0(xi2) + c11*phi2(xi1)*phi1(xi2) + c12*phi3(xi1)*phi0(xi2) + c13*phi3(xi1)*phi1(xi2))*x1 \
    + (c20*phi0(xi1)*phi2(xi2) + c21*phi0(xi1)*phi3(xi2) + c22*phi1(xi1)*phi2(xi2) + c23*phi1(xi1)*phi3(xi2))*x2 \
    + (c30*phi2(xi1)*phi2(xi2) + c31*phi2(xi1)*phi3(xi2) + c32*phi3(xi1)*phi2(xi2) + c33*phi3(xi1)*phi3(xi2))*x3

print ""
print "x(t) =",factor(x_2(xi1,xi2),[xi1,xi2])
p1 = factor(x_1(xi1,xi2), [xi1,xi2])
p2 = factor(x_2(xi1,xi2), [xi1,xi2])
print ""
#print p1
print ""

print "factors:"
print [xi1**ex1*xi2**ex2 for ex1 in range(0,4) for ex2 in range(0,4)]
coeff1 = [p1.coeff(xi1**ex1*xi2**ex2) for ex1 in range(0,4) for ex2 in range(0,4)]
coeff2 = [p2.coeff(xi1**ex1*xi2**ex2) for ex1 in range(0,4) for ex2 in range(0,4)]

print "coefficients of x_1: ", coeff1
print "coefficients of x_2: ", coeff2

print "solve for c's"
result = solve([coeff1[i] - coeff2[i] for i in range(16)],[c00,c01,c02,c03,c10,c11,c12,c13,c20,c21,c22,c23,c30,c31,c32,c33])
print result
print "test:"
p2 = x_2(xi1,xi2)
print ""
print factor(p2.subs(result),[xi1,xi2])
print ""
print "conclusion: to represent a 1D linear geometry in an element between xa and xb, use the following nodal dofs:"
print "node0: dof0 (c0):",result[c00]
print "node0: dof1 (c1):",result[c01]
print "node1: dof0 (c2):",result[c02]
print "node1: dof1 (c3):",result[c03]

###########################################################################
### 3D element contributions for linear Lagrange ansatz function (phi)  ###

def phi0(xi1,xi2,xi3):
  return (1-xi1)*(1-xi2)*(1-xi3)
  
def phi1(xi1,xi2,xi3):
  return xi1*(1-xi2)*(1-xi3)

def phi2(xi1,xi2,xi3):
  return (1-xi1)*xi2*(1-xi3)

def phi3(xi1,xi2,xi3):
  return xi1*xi2*(1-xi3)
  
def phi4(xi1,xi2,xi3):
  return (1-xi1)*(1-xi2)*xi3
  
def phi5(xi1,xi2,xi3):
  return xi1*(1-xi2)*xi3

def phi6(xi1,xi2,xi3):
  return (1-xi1)*xi2*xi3

def phi7(xi1,xi2,xi3):
  return xi1*xi2*xi3

print ""
print "3D stencil"
print "=========="

print diff(phi0(xi1,xi2,xi3),xi1).subs(xi1,0.5).subs(xi2,0.5).subs(xi3,0.5)
print diff(phi0(xi1,xi2,xi3),xi2).subs(xi1,0.5).subs(xi2,0.5).subs(xi3,0.5)
print diff(phi0(xi1,xi2,xi3),xi3).subs(xi1,0.5).subs(xi2,0.5).subs(xi3,0.5)

a1 = diff(phi0(xi1,xi2,xi3),xi1)
a2 = diff(phi0(xi1,xi2,xi3),xi2)
a3 = diff(phi0(xi1,xi2,xi3),xi3)

integrand = a1*a1 + a2*a2 + a3*a3
print "integrand: ",integrand
print "integrand: ",integrand.subs(xi1,0.5).subs(xi2,0.5).subs(xi3,0.5)
print "symbolic ", integrate(integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))


# bottom lower left node
integrand = diff(phi0(xi1,xi2,xi3),xi1)*diff(phi0(xi1,xi2,xi3),xi1) + diff(phi0(xi1,xi2,xi3),xi2)*diff(phi0(xi1,xi2,xi3),xi2) + diff(phi0(xi1,xi2,xi3),xi3)*diff(phi0(xi1,xi2,xi3),xi3)
entry000 = integrate(-integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # -1/3

# bottom lower right node
integrand = diff(phi0(xi1,xi2,xi3),xi1)*diff(phi1(xi1,xi2,xi3),xi1) + diff(phi0(xi1,xi2,xi3),xi2)*diff(phi1(xi1,xi2,xi3),xi2) + diff(phi0(xi1,xi2,xi3),xi3)*diff(phi1(xi1,xi2,xi3),xi3)
entry001 = integrate(-integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 0

# bottom upper left node
integrand = diff(phi0(xi1,xi2,xi3),xi1)*diff(phi2(xi1,xi2,xi3),xi1) + diff(phi0(xi1,xi2,xi3),xi2)*diff(phi2(xi1,xi2,xi3),xi2) + diff(phi0(xi1,xi2,xi3),xi3)*diff(phi2(xi1,xi2,xi3),xi3)
entry010 = integrate(-integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 0

# bottom upper right node
integrand = diff(phi0(xi1,xi2,xi3),xi1)*diff(phi3(xi1,xi2,xi3),xi1) + diff(phi0(xi1,xi2,xi3),xi2)*diff(phi3(xi1,xi2,xi3),xi2) + diff(phi0(xi1,xi2,xi3),xi3)*diff(phi3(xi1,xi2,xi3),xi3)
entry011 = integrate(-integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 1/12

# top lower left node
integrand = diff(phi0(xi1,xi2,xi3),xi1)*diff(phi4(xi1,xi2,xi3),xi1) + diff(phi0(xi1,xi2,xi3),xi2)*diff(phi4(xi1,xi2,xi3),xi2) + diff(phi0(xi1,xi2,xi3),xi3)*diff(phi4(xi1,xi2,xi3),xi3)
entry100 = integrate(-integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 0

# top lower right node
integrand = diff(phi0(xi1,xi2,xi3),xi1)*diff(phi5(xi1,xi2,xi3),xi1) + diff(phi0(xi1,xi2,xi3),xi2)*diff(phi5(xi1,xi2,xi3),xi2) + diff(phi0(xi1,xi2,xi3),xi3)*diff(phi5(xi1,xi2,xi3),xi3)
entry101 = integrate(-integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 1/12

# top upper left node
integrand = diff(phi0(xi1,xi2,xi3),xi1)*diff(phi6(xi1,xi2,xi3),xi1) + diff(phi0(xi1,xi2,xi3),xi2)*diff(phi6(xi1,xi2,xi3),xi2) + diff(phi0(xi1,xi2,xi3),xi3)*diff(phi6(xi1,xi2,xi3),xi3)
entry110 = integrate(-integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 1/12

# top upper right node
integrand = diff(phi0(xi1,xi2,xi3),xi1)*diff(phi7(xi1,xi2,xi3),xi1) + diff(phi0(xi1,xi2,xi3),xi2)*diff(phi7(xi1,xi2,xi3),xi2) + diff(phi0(xi1,xi2,xi3),xi3)*diff(phi7(xi1,xi2,xi3),xi3)
entry111 = integrate(-integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 1/12

print "element contribution:"
print "center: [",entry010,entry011,"]"
print "        [",entry000,entry001,"]"
print ""
print "top:    [",entry110,entry111,"]"
print "        [",entry100,entry101,"]"
print ""
print "nodal stencil:"
print "bottom: [",entry111,2*entry110,entry111,"]"
print "        [",2*entry101,4*entry100,2*entry101,"]"
print "        [",entry111,2*entry110,entry111,"]"
print ""
print "center: [",2*entry011,4*entry010,2*entry011,"]"
print "        [",4*entry001,8*entry000,4*entry001,"]"
print "        [",2*entry011,4*entry010,2*entry011,"]"
print ""
print "top:    [",entry111,2*entry110,entry111,"]"
print "        [",2*entry101,4*entry100,2*entry101,"]"
print "        [",entry111,2*entry110,entry111,"]"

# the element contribution stencil is
# center: [   0     1/12 ]
#         [ _-1/3_  0    ]
#
# top:    [   1/12  1/12 ]
#         [   0     1/12 ]

# the full stencil is therefore
# bottom: [   1/3    1/6   1/3  ]         [   1   2   1  ]   
#         [   1/6    0     1/6  ]         [   2   0   2  ]   
#         [   1/3    1/6   1/3  ]         [   1   2   1  ]   
#                                                                  
# center: [   1/6    0     1/6  ]         [   2   0   2  ]   
#         [   0    _-8/3_  0    ] =  1/12*[   0 _-32_ 0  ]
#         [   1/6    0     1/6  ]         [   2   0   2  ]   
#                                                                  
# top:    [   1/3    1/6   1/3  ]         [   1   2   1  ]   
#         [   1/6    0     1/6  ]         [   2   0   2  ]   
#         [   1/3    1/6   1/3  ]         [   1   2   1  ]   


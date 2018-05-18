#!/usr/bin/pxi2thon
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *

# compute stencil for FE rhs (int phi_i phi_j)
#
xi,xi1,xi2,xi3 = symbols('xi xi1 xi2 xi3')

### 1D element contributions for linear Lagrange ansatz function (phi)

def phi0(xi):
  return (1-xi)
  
def phi1(xi):
  return xi
  
print "1D stencil"
print "=========="

# left node
integrand = phi0(xi)**2
entry0 = integrate(integrand,(xi,0,1))

# right node
integrand = phi0(xi)*phi1(xi)
entry1 = integrate(integrand,(xi,0,1))

print "element contribution:"
print "[",entry0,entry1,"]"

print "nodal stencil:"
print "[",entry1,2*entry0,entry1,"]"

### 2D element contributions for linear Lagrange ansatz function (phi)

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
integrand = phi0(xi1,xi2)**2
entry00 = integrate(integrand,(xi1,0,1),(xi2,0,1))      # -2/3

# lower right node
integrand = phi0(xi1,xi2)*phi1(xi1,xi2)
entry01 = integrate(integrand,(xi1,0,1),(xi2,0,1))      # 1/6

# upper left node
integrand = phi0(xi1,xi2)*phi2(xi1,xi2)
entry10 = integrate(integrand,(xi1,0,1),(xi2,0,1))      # 1/6

# upper right node
integrand = phi0(xi1,xi2)*phi3(xi1,xi2)
entry11 = integrate(integrand,(xi1,0,1),(xi2,0,1))      # 1/3

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


### 3D element contributions for linear Lagrange ansatz function (phi)

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

# bottom lower left node
integrand = phi0(xi1,xi2,xi3)*phi0(xi1,xi2,xi3)
entry000 = integrate(integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # -1/3

# bottom lower right node
integrand = phi0(xi1,xi2,xi3)*phi1(xi1,xi2,xi3)
entry001 = integrate(integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 0

# bottom upper left node
integrand = phi0(xi1,xi2,xi3)*phi2(xi1,xi2,xi3)
entry010 = integrate(integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 0

# bottom upper right node
integrand = phi0(xi1,xi2,xi3)*phi3(xi1,xi2,xi3)
entry011 = integrate(integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 1/12

# top lower left node
integrand = phi0(xi1,xi2,xi3)*phi4(xi1,xi2,xi3)
entry100 = integrate(integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 0

# top lower right node
integrand = phi0(xi1,xi2,xi3)*phi5(xi1,xi2,xi3)
entry101 = integrate(integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 1/12

# top upper left node
integrand = phi0(xi1,xi2,xi3)*phi6(xi1,xi2,xi3)
entry110 = integrate(integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 1/12

# top upper right node
integrand = phi0(xi1,xi2,xi3)*phi7(xi1,xi2,xi3)
entry111 = integrate(integrand,(xi1,0,1),(xi2,0,1),(xi3,0,1))      # 1/12

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


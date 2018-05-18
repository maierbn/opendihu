#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
from sympy.matrices import *
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False, no_global=True)
from sympy.printing import cxxcode

if False:
  m11,m12,m13 = symbols('m11,m12,m13')
  m21,m22,m23 = symbols('m21,m22,m23')
  m31,m32,m33 = symbols('m31,m32,m33')

  # compute inverse of M
  M=Matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])
  print "matrix:"
  print M
  print ""
  print "inv:",simplify(M.inv())
  # inv = 1/det * adj

  print "adj:",simplify(M.adjugate())
   
  print "det:",simplify(M.det())
  print "T: 1/(det**2)*",simplify(M.adjugate()*M.adjugate().transpose())
   
  # compute inverse of a symmetric matrix
  M=Matrix([[m11,m21,m31],[m21,m22,m32],[m31,m32,m33]])
  print ""
  print "symmetric matrix:"
  print M
  print ""
  print "inv:",simplify(M.inv())
  print ""
  print "adj:",simplify(M.adjugate())
  print ""
  print "det:",simplify(M.det())
  print ""
  print ""
   
  # 3D
  # compute T
  l1,l2 = symbols('l1,l2')
  M2 = Matrix([[l1, -l1*cos(a)/sin(a)],[0,l2/sin(a)]])
  M2inv = M2.inv()
  print "M2:",M2
  print "inv:",M2inv
  print "invT:",M2inv.transpose()
  print "T:",M2inv*M2inv.transpose()

# mapping from parameter space to world space

xi1,xi2,xi3 = symbols('xi1,xi2,xi3', positive=True, real=True)
xp1,x11,x21,x31,x41,x51,x61,x71,x81 = symbols('xp1,x11,x21,x31,x41,x51,x61,x71,x81', real=True)
xp2,x12,x22,x32,x42,x52,x62,x72,x82 = symbols('xp2,x12,x22,x32,x42,x52,x62,x72,x82', real=True)
xp3,x13,x23,x33,x43,x53,x63,x73,x83 = symbols('xp3,x13,x23,x33,x43,x53,x63,x73,x83', real=True)

# 1D linear
Phi1 = (1-xi1)*x11 + xi1*x21
Phi2 = (1-xi1)*x12 + xi1*x22
Phi3 = (1-xi1)*x13 + xi1*x22

print "1D linear"
print solve(Phi1 - xp1, xi1)

# 2D linear
Phi1 = (1-xi1)*(1-xi2)*x11 + xi1*(1-xi2)*x21 + (1-xi1)*xi2*x31 + xi1*xi2*x41
Phi2 = (1-xi1)*(1-xi2)*x12 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x32 + xi1*xi2*x42
Phi3 = (1-xi1)*(1-xi2)*x13 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x33 + xi1*xi2*x43

print "2D linear"
solution = solve([Phi1 - xp1, Phi2 - xp2], [xi1,xi2])
print len(solution)," solution(s)"
solution = simplify(solution)
print solution
print ""

expr = collect(solution[0][0], [xp1,xp2])
print "xi1=",expr
print(cxxcode(expr, standard='C++11'))
print ""

expr = collect(solution[0][1], [xp1,xp2])
print "xi2=",expr
print(cxxcode(expr, standard='C++11'))
print ""

expr = simplify(solve(Phi1 - xp1, xi2))[0]
print "compute value of xi1 and then compute xi2=",expr
print(cxxcode(expr, standard='C++11'))

# 3D linear
Phi1 = (1-xi3)*((1-xi1)*(1-xi2)*x11 + xi1*(1-xi2)*x21 + (1-xi1)*xi2*x31 + xi1*xi2*x41) + xi3*((1-xi1)*(1-xi2)*x51 + xi1*(1-xi2)*x61 + (1-xi1)*xi2*x71 + xi1*xi2*x81)
Phi2 = (1-xi3)*((1-xi1)*(1-xi2)*x12 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x32 + xi1*xi2*x42) + xi3*((1-xi1)*(1-xi2)*x52 + xi1*(1-xi2)*x62 + (1-xi1)*xi2*x72 + xi1*xi2*x82)
Phi3 = (1-xi3)*((1-xi1)*(1-xi2)*x13 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x33 + xi1*xi2*x43) + xi3*((1-xi1)*(1-xi2)*x53 + xi1*(1-xi2)*x62 + (1-xi1)*xi2*x73 + xi1*xi2*x83)

print "3D linear"
solution = solve([Phi1 - xp1, Phi2 - xp2, Phi3 - xp3], [xi1,xi2,xi3])
print len(solution)," solution(s)"
print solution
print "----------"
solution = simplify(solution)
print solution
print ""

expr = collect(solution[0][0], [xp1,xp2,xp3])
print "xi1=",expr
print(cxxcode(expr, standard='C++11'))
print ""

expr = collect(solution[0][1], [xp1,xp2,xp3])
print "xi2=",expr
print(cxxcode(expr, standard='C++11'))
print ""

expr = collect(solution[0][2], [xp1,xp2,xp3])
print "xi3=",expr
print(cxxcode(expr, standard='C++11'))
print ""

expr = collect(solve([Phi1 - xp1, Phi2 - xp2], [xi2,xi3]), [xp1,xp2,xp3])[0]
print "compute value of xi1 and then compute xi2=",expr
print(cxxcode(expr, standard='C++11'))
print ""

expr = collect(solve(Phi1 - xp1, xi3), [xp1,xp2,xp3])
print "then compute xi3=",expr[0]
print(cxxcode(expr, standard='C++11'))

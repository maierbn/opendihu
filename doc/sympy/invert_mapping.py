#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
from sympy.matrices import *
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False, no_global=True)
from sympy.printing import cxxcode
import sys

if True:
  m11,m12,m13 = symbols('m11,m12,m13')
  m21,m22,m23 = symbols('m21,m22,m23')
  m31,m32,m33 = symbols('m31,m32,m33')

  # compute inverse of M
  M=Matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])
  print "matrix:"
  print M
  print ""
  print M[0],M[3],M[6]
  print ""
  
  
  print "inv:",simplify(M.inv())
  # inv = 1/det * adj

  print "adj:",simplify(M.adjugate())
   
  print "det:",simplify(M.det())
  print "T: 1/(det**2)*",simplify(M.adjugate()*M.adjugate().transpose())
   
  sys.exit(0)
   
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
if False:
  Phi1 = (1-xi1)*x11 + xi1*x21
  Phi2 = (1-xi1)*x12 + xi1*x22
  Phi3 = (1-xi1)*x13 + xi1*x22

  print "1D linear"
  print solve(Phi1 - xp1, xi1)

# 2D linear quadrilateral
if True:
  Phi1 = (1-xi1)*(1-xi2)*x11 + xi1*(1-xi2)*x21 + (1-xi1)*xi2*x31 + xi1*xi2*x41
  Phi2 = (1-xi1)*(1-xi2)*x12 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x32 + xi1*xi2*x42
  Phi3 = (1-xi1)*(1-xi2)*x13 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x33 + xi1*xi2*x43

  print "2D linear"
  solution = solve([Phi1 - xp1, Phi2 - xp2], [xi1,xi2])
  print len(solution)," solution(s)"
  solution = simplify(solution)
  print solution
  print ""

  # the 2nd solution is the correct one
  expr = collect(solution[1][0], [xp1,xp2])
  print "xi1=",expr
  print(cxxcode(expr, standard='C++11'))
  print ""

  expr = collect(solution[1][1], [xp1,xp2])
  print "xi2=",expr
  #print(cxxcode(expr, standard='C++11'))
  print ""

  expr = simplify(solve(Phi1 - xp1, xi2))[0]
  print "compute value of xi1 and then compute xi2=",expr
  print(cxxcode(expr, standard='C++11'))

#2D linear simplex
if False:
  Phi1 = (1-xi1-xi2)*x11 + xi1*x21 + xi2*x31
  Phi2 = (1-xi1-xi2)*x12 + xi1*x22 + xi2*x32
  Phi3 = (1-xi1-xi2)*x13 + xi1*x23 + xi2*x33
  
  print "2D simplex"
  solution = solve([Phi1 - xp1, Phi2 - xp2], [xi1,xi2])
  solution = simplify(solution)
  print solution
  print ""
  
  expr = collect(solution[xi1], [xp1,xp2])
  print "xi1=",expr
  print "c++:", cxxcode(expr, standard='C++11')
  print ""

  expr = collect(solution[xi2], [xp1,xp2])
  print "xi2=",expr
  print "c++:", cxxcode(expr, standard='C++11')
  print ""

  expr = collect(solve([Phi1 - xp1], [xi2])[xi2], [xp1,xp2])
  print "compute value of xi1 and then compute xi2=",expr
  print "c++:", cxxcode(expr, standard='C++11')
  print ""
  
#3D linear simplex
if False:
  Phi1 = (1-xi1-xi2-xi3)*x11 + xi1*x21 + xi2*x31 + xi3*x41
  Phi2 = (1-xi1-xi2-xi3)*x12 + xi1*x22 + xi2*x32 + xi3*x42
  Phi3 = (1-xi1-xi2-xi3)*x13 + xi1*x23 + xi2*x33 + xi3*x43
  
  print "3D simplex"
  solution = solve([Phi1 - xp1, Phi2 - xp2, Phi3 - xp3], [xi1,xi2,xi3])
  solution = simplify(solution)
  print solution
  print ""
  
  expr = collect(solution[xi1], [xp1,xp2,xp3])
  print "xi1=",expr
  print "c++:", cxxcode(expr, standard='C++11')
  print ""

  expr = collect(solution[xi2], [xp1,xp2,xp3])
  print "xi2=",expr
  print "c++:", cxxcode(expr, standard='C++11')
  print ""
  
  expr = collect(solution[xi3], [xp1,xp2,xp3])
  print "xi3=",expr
  print "c++:", cxxcode(expr, standard='C++11')
  print ""

  expr = collect(solve([Phi1 - xp1, Phi2 - xp2], [xi2,xi3])[xi2], [xp1,xp2,xp3])
  print "compute value of xi1 and then compute xi2=",expr
  print "c++:", cxxcode(expr, standard='C++11')
  print ""
  
  expr = collect(solve([Phi1 - xp1], [xi3])[xi3], [xp1,xp2,xp3])
  print "then compute xi3=",expr
  print "c++:", cxxcode(expr, standard='C++11')
  print ""

if False:
  #it doesn't work like this, kl
  print "3D simplex, the inner, big one when decomposing the hexaeder into 5 tetrahedrons"
  
  t1,t2,t3,t4,t5 = symbols('t1,t2,t3,t4,t5', real=True)
  xia1,xia4,xia7 = symbols('xia1,xia4,xia7', positive=True, real=True)
  
  xp1,xp2,xp3 = symbols('xp1,xp2,xp3', real=True)
  
  #xp1 = (1-xia1-xia4-xia7)*x21 + xia1*x11 + xia4*x41 + xia7*x71
  #xp2 = (1-xia1-xia4-xia7)*x22 + xia1*x12 + xia4*x42 + xia7*x72
  #xp3 = (1-xia1-xia4-xia7)*x23 + xia1*x13 + xia4*x43 + xia7*x73

  Phi1 = (1-xi3)*((1-xi1)*(1-xi2)*x11 + xi1*(1-xi2)*x21 + (1-xi1)*xi2*x31 + xi1*xi2*x41) + xi3*((1-xi1)*(1-xi2)*x51 + xi1*(1-xi2)*x61 + (1-xi1)*xi2*x71 + xi1*xi2*x81)
  Phi2 = (1-xi3)*((1-xi1)*(1-xi2)*x12 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x32 + xi1*xi2*x42) + xi3*((1-xi1)*(1-xi2)*x52 + xi1*(1-xi2)*x62 + (1-xi1)*xi2*x72 + xi1*xi2*x82)
  Phi3 = (1-xi3)*((1-xi1)*(1-xi2)*x13 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x33 + xi1*xi2*x43) + xi3*((1-xi1)*(1-xi2)*x53 + xi1*(1-xi2)*x62 + (1-xi1)*xi2*x73 + xi1*xi2*x83)
 
  #xi1sol = solve([xp1-Phi1], [xi1])[xi1]
  #print "xi1: ", xi1sol
  #print ""
  
  #Phi2sol = simplify(Phi2.subs(xi1,xi1sol))
  #print "Phi2sol:", Phi2sol
  #print ""
  #Phi2sol = factor(Phi2sol)
  #print "Phi2sol: ",Phi2sol
  #print ""
  
  #Phi2sol = collect(Phi2sol, xi2)
  #print "Phi2sol: ",Phi2sol
  #print ""
  
  print "sustitute Phi2sol"
  Phi2sol = (t1 + xi2*t2 + xi2**2*t3) / (t4 + xi2*t5)
  
  t1 = -x11*x12*xi3*xia1 + x11*x12*xia1 + x11*x22*xi3**2 + x11*x22*xi3*xia1 - 2*x11*x22*xi3 - x11*x22*xia1 + x11*x22 + x11*x52*xi3*xia1 - x11*x62*xi3**2 - x11*x62*xi3*xia1 + x11*x62*xi3 - x12*x21*xi3**2 + x12*x21*xi3*xia1 + x12*x21*xi3*xia4 + x12*x21*xi3*xia7 + x12*x21*xi3 - x12*x21*xia1 - x12*x21*xia4 - x12*x21*xia7 - x12*x41*xi3*xia4 + x12*x41*xia4 + x12*x61*xi3**2 - x12*x61*xi3 - x12*x71*xi3*xia7 + x12*x71*xia7 - x21*x22*xi3*xia1 - x21*x22*xi3*xia4 - x21*x22*xi3*xia7 + x21*x22*xi3 + x21*x22*xia1 + x21*x22*xia4 + x21*x22*xia7 - x21*x22 + x21*x52*xi3**2 - x21*x52*xi3*xia1 - x21*x52*xi3*xia4 - x21*x52*xi3*xia7 + x21*x62*xi3*xia1 + x21*x62*xi3*xia4 + x21*x62*xi3*xia7 - x21*x62*xi3 + x22*x41*xi3*xia4 - x22*x41*xia4 - x22*x51*xi3**2 + x22*x51*xi3 + x22*x71*xi3*xia7 - x22*x71*xia7 + x41*x52*xi3*xia4 - x41*x62*xi3*xia4 + x51*x62*xi3**2 - x52*x61*xi3**2 + x52*x71*xi3*xia7 - x62*x71*xi3*xia7 
  t2 = x11*x12*xi3*xia1 - x11*x12*xia1 - 2*x11*x22*xi3**2 - x11*x22*xi3*xia1 + 4*x11*x22*xi3 + x11*x22*xia1 - 2*x11*x22 - x11*x32*xi3*xia1 + x11*x32*xia1 + x11*x42*xi3**2 + x11*x42*xi3*xia1 - 2*x11*x42*xi3 - x11*x42*xia1 + x11*x42 - x11*x52*xi3*xia1 + 2*x11*x62*xi3**2 + x11*x62*xi3*xia1 - 2*x11*x62*xi3 + x11*x72*xi3*xia1 - x11*x82*xi3**2 - x11*x82*xi3*xia1 + x11*x82*xi3 + 2*x12*x21*xi3**2 - x12*x21*xi3*xia1 - x12*x21*xi3*xia4 - x12*x21*xi3*xia7 - 3*x12*x21*xi3 + x12*x21*xia1 + x12*x21*xia4 + x12*x21*xia7 + x12*x21 - x12*x41*xi3**2 + x12*x41*xi3*xia4 + 2*x12*x41*xi3 - x12*x41*xia4 - x12*x41 - 2*x12*x61*xi3**2 + 2*x12*x61*xi3 + x12*x71*xi3*xia7 - x12*x71*xia7 + x12*x81*xi3**2 - x12*x81*xi3 + x21*x22*xi3*xia1 + x21*x22*xi3*xia4 + x21*x22*xi3*xia7 - x21*x22*xi3 - x21*x22*xia1 - x21*x22*xia4 - x21*x22*xia7 + x21*x22 - x21*x32*xi3**2 + x21*x32*xi3*xia1 + x21*x32*xi3*xia4 + x21*x32*xi3*xia7 + x21*x32*xi3 - x21*x32*xia1 - x21*x32*xia4 - x21*x32*xia7 - x21*x42*xi3*xia1 - x21*x42*xi3*xia4 - x21*x42*xi3*xia7 + x21*x42*xi3 + x21*x42*xia1 + x21*x42*xia4 + x21*x42*xia7 - x21*x42 - 2*x21*x52*xi3**2 + x21*x52*xi3*xia1 + x21*x52*xi3*xia4 + x21*x52*xi3*xia7 + x21*x52*xi3 - x21*x62*xi3*xia1 - x21*x62*xi3*xia4 - x21*x62*xi3*xia7 + x21*x62*xi3 + x21*x72*xi3**2 - x21*x72*xi3*xia1 - x21*x72*xi3*xia4 - x21*x72*xi3*xia7 + x21*x82*xi3*xia1 + x21*x82*xi3*xia4 + x21*x82*xi3*xia7 - x21*x82*xi3 + x22*x31*xi3**2 - 2*x22*x31*xi3 + x22*x31 - x22*x41*xi3*xia4 + x22*x41*xia4 + 2*x22*x51*xi3**2 - 2*x22*x51*xi3 - x22*x71*xi3**2 - x22*x71*xi3*xia7 + x22*x71*xi3 + x22*x71*xia7 - x31*x62*xi3**2 + x31*x62*xi3 - x32*x41*xi3*xia4 + x32*x41*xia4 + x32*x61*xi3**2 - x32*x61*xi3 - x32*x71*xi3*xia7 + x32*x71*xia7 + x41*x42*xi3*xia4 - x41*x42*xia4 + x41*x52*xi3**2 - x41*x52*xi3*xia4 - x41*x52*xi3 + x41*x62*xi3*xia4 + x41*x72*xi3*xia4 - x41*x82*xi3*xia4 - x42*x51*xi3**2 + x42*x51*xi3 + x42*x71*xi3*xia7 - x42*x71*xia7 - 2*x51*x62*xi3**2 + x51*x82*xi3**2 + 2*x52*x61*xi3**2 - x52*x71*xi3*xia7 - x52*x81*xi3**2 - x61*x72*xi3**2 + x62*x71*xi3**2 + x62*x71*xi3*xia7 + x71*x72*xi3*xia7 - x71*x82*xi3*xia7
  t3 = x11*x22*xi3**2 - 2*x11*x22*xi3 + x11*x22 - x11*x42*xi3**2 + 2*x11*x42*xi3 - x11*x42 - x11*x62*xi3**2 + x11*x62*xi3 + x11*x82*xi3**2 - x11*x82*xi3 - x12*x21*xi3**2 + 2*x12*x21*xi3 - x12*x21 + x12*x41*xi3**2 - 2*x12*x41*xi3 + x12*x41 + x12*x61*xi3**2 - x12*x61*xi3 - x12*x81*xi3**2 + x12*x81*xi3 + x21*x32*xi3**2 - 2*x21*x32*xi3 + x21*x32 + x21*x52*xi3**2 - x21*x52*xi3 - x21*x72*xi3**2 + x21*x72*xi3 - x22*x31*xi3**2 + 2*x22*x31*xi3 - x22*x31 - x22*x51*xi3**2 + x22*x51*xi3 + x22*x71*xi3**2 - x22*x71*xi3 + x31*x42*xi3**2 - 2*x31*x42*xi3 + x31*x42 + x31*x62*xi3**2 - x31*x62*xi3 - x31*x82*xi3**2 + x31*x82*xi3 - x32*x41*xi3**2 + 2*x32*x41*xi3 - x32*x41 - x32*x61*xi3**2 + x32*x61*xi3 + x32*x81*xi3**2 - x32*x81*xi3 - x41*x52*xi3**2 + x41*x52*xi3 + x41*x72*xi3**2 - x41*x72*xi3 + x42*x51*xi3**2 - x42*x51*xi3 - x42*x71*xi3**2 + x42*x71*xi3 + x51*x62*xi3**2 - x51*x82*xi3**2 - x52*x61*xi3**2 + x52*x81*xi3**2 + x61*x72*xi3**2 - x62*x71*xi3**2 + x71*x82*xi3**2 - x72*x81*xi3**2
  t4 = -x11*xi3 + x11 + x21*xi3 - x21 + x51*xi3 - x61*xi3
  t5 = x11*xi3 - x11 - x21*xi3 + x21 - x31*xi3 + x31 + x41*xi3 - x41 - x51*xi3 + x61*xi3 + x71*xi3 - x81*xi3
  
  print solve([xp2-Phi2sol], [xi2])[0]
  
  xi2_a = (-t2 + t5*xp2 - sqrt(-4*t1*t3 + t2**2 - 2*t2*t5*xp2 + 4*t3*t4*xp2 + t5**2*xp2**2))/(2*t3)
  xi2_b = (-t2 + t5*xp2 + sqrt(-4*t1*t3 + t2**2 - 2*t2*t5*xp2 + 4*t3*t4*xp2 + t5**2*xp2**2))/(2*t3)
  
  
  
  
  xi2 = solve([xp2-Phi2sol], [xi2])[xi2]
  print "xi2: ", xi2
  print ""
  
  xi3 = solve([xp3-Phi3], [xi3])[xi3]
  print "xi3: ",xi3
  print ""
  
 
  print solve([xp1-Phi1, xp2-Phi2, xp3-Phi3], [xi1,xi2,xi3])

if True:
  
  T = Matrix([[x11-x41,x21-x41,x31-x41],[x12-x42,x22-x42,x32-x42],[x13-x43,x23-x43,x33-x43]])
  tadj = T.adjugate()
  inv = T.inv()
  print "adj:",tadj
  print "inv:",inv
  
  xi1 = tadj[0]*(xp1-x41) + tadj[1]*(xp2-x42) + tadj[2]*(xp3-x43)
  xi2 = tadj[3]*(xp1-x41) + tadj[4]*(xp2-x42) + tadj[5]*(xp3-x43)
  xi3 = tadj[6]*(xp1-x41) + tadj[7]*(xp2-x42) + tadj[8]*(xp3-x43)
  
  print "det =",T.det()
  print "xi1 = 1/det *",xi1
  print "xi2 = 1/det *",xi2
  print "xi3 = 1/det *",xi3

# 3D linear
if False:
  
  Phi1 = (1-xi3)*((1-xi1)*(1-xi2)*x11 + xi1*(1-xi2)*x21 + (1-xi1)*xi2*x31 + xi1*xi2*x41) + xi3*((1-xi1)*(1-xi2)*x51 + xi1*(1-xi2)*x61 + (1-xi1)*xi2*x71 + xi1*xi2*x81)
  
  print "solve Phi1 for xi2:"
  xi2sol = solve(Phi1 - xp1,xi2)
  print "xi2(xi1,xi3)=",xi2sol
  print ""
  
  Phi1 = (1-xi1)*(1-xi2)*x11 + xi1*(1-xi2)*x21 + (1-xi1)*xi2*x31 + xi1*xi2*x41
  Phi2 = (1-xi1)*(1-xi2)*x12 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x32 + xi1*xi2*x42

  print "(2D linear)"
  #solution12 = solve([Phi1 - xp1, Phi2 - xp2], [xi1,xi2])
  #print len(solution12)," solution(s)"
  #solution12 = simplify(solution12)
  #print solution12
  print ""

  # the second solution is the correct one
  #expr = solution12[1][0]
  #expr = (2*x11*x32 - x11*x42 - x11*xp2 - 2*x12*x31 + x12*x41 + x12*xp1 - x21*x32 + x21*xp2 + x22*x31 - x22*xp1 + x31*xp2 - x32*xp1 - x41*xp2 + x42*xp1 + sqrt(x11**2*x42**2 - 2*x11**2*x42*xp2 + x11**2*xp2**2 - 2*x11*x12*x41*x42 + 2*x11*x12*x41*xp2 + 2*x11*x12*x42*xp1 - 2*x11*x12*xp1*xp2 - 2*x11*x21*x32*x42 + 2*x11*x21*x32*xp2 + 2*x11*x21*x42*xp2 - 2*x11*x21*xp2**2 - 2*x11*x22*x31*x42 + 2*x11*x22*x31*xp2 + 4*x11*x22*x32*x41 - 4*x11*x22*x32*xp1 - 4*x11*x22*x41*xp2 + 2*x11*x22*x42*xp1 + 2*x11*x22*xp1*xp2 + 2*x11*x31*x42*xp2 - 2*x11*x31*xp2**2 - 4*x11*x32*x41*xp2 + 2*x11*x32*x42*xp1 + 2*x11*x32*xp1*xp2 + 2*x11*x41*x42*xp2 + 2*x11*x41*xp2**2 - 2*x11*x42**2*xp1 - 2*x11*x42*xp1*xp2 + x12**2*x41**2 - 2*x12**2*x41*xp1 + x12**2*xp1**2 + 4*x12*x21*x31*x42 - 4*x12*x21*x31*xp2 - 2*x12*x21*x32*x41 + 2*x12*x21*x32*xp1 + 2*x12*x21*x41*xp2 - 4*x12*x21*x42*xp1 + 2*x12*x21*xp1*xp2 - 2*x12*x22*x31*x41 + 2*x12*x22*x31*xp1 + 2*x12*x22*x41*xp1 - 2*x12*x22*xp1**2 + 2*x12*x31*x41*xp2 - 4*x12*x31*x42*xp1 + 2*x12*x31*xp1*xp2 + 2*x12*x32*x41*xp1 - 2*x12*x32*xp1**2 - 2*x12*x41**2*xp2 + 2*x12*x41*x42*xp1 - 2*x12*x41*xp1*xp2 + 2*x12*x42*xp1**2 + x21**2*x32**2 - 2*x21**2*x32*xp2 + x21**2*xp2**2 - 2*x21*x22*x31*x32 + 2*x21*x22*x31*xp2 + 2*x21*x22*x32*xp1 - 2*x21*x22*xp1*xp2 + 2*x21*x31*x32*xp2 - 4*x21*x31*x42*xp2 + 2*x21*x31*xp2**2 - 2*x21*x32**2*xp1 + 2*x21*x32*x41*xp2 + 2*x21*x32*x42*xp1 - 2*x21*x32*xp1*xp2 - 2*x21*x41*xp2**2 + 2*x21*x42*xp1*xp2 + x22**2*x31**2 - 2*x22**2*x31*xp1 + x22**2*xp1**2 - 2*x22*x31**2*xp2 + 2*x22*x31*x32*xp1 + 2*x22*x31*x41*xp2 + 2*x22*x31*x42*xp1 - 2*x22*x31*xp1*xp2 - 4*x22*x32*x41*xp1 + 2*x22*x32*xp1**2 + 2*x22*x41*xp1*xp2 - 2*x22*x42*xp1**2 + x31**2*xp2**2 - 2*x31*x32*xp1*xp2 - 2*x31*x41*xp2**2 + 2*x31*x42*xp1*xp2 + x32**2*xp1**2 + 2*x32*x41*xp1*xp2 - 2*x32*x42*xp1**2 + x41**2*xp2**2 - 2*x41*x42*xp1*xp2 + x42**2*xp1**2))/(2*(x11*x32 - x11*x42 - x12*x31 + x12*x41 - x21*x32 + x21*x42 + x22*x31 - x22*x41))  
  
  # xi1
  expr = (2*x11*x32 - x11*x42 - 2*x12*x31 + x12*x41 - x21*x32 + x22*x31 + xp1*(x12 - x22 - x32 + x42) + xp2*(-x11 + x21 + x31 - x41) + sqrt(x11**2*x42**2 - 2*x11**2*x42*xp2 + x11**2*xp2**2 - 2*x11*x12*x41*x42 + 2*x11*x12*x41*xp2 + 2*x11*x12*x42*xp1 - 2*x11*x12*xp1*xp2 - 2*x11*x21*x32*x42 + 2*x11*x21*x32*xp2 + 2*x11*x21*x42*xp2 - 2*x11*x21*xp2**2 - 2*x11*x22*x31*x42 + 2*x11*x22*x31*xp2 + 4*x11*x22*x32*x41 - 4*x11*x22*x32*xp1 - 4*x11*x22*x41*xp2 + 2*x11*x22*x42*xp1 + 2*x11*x22*xp1*xp2 + 2*x11*x31*x42*xp2 - 2*x11*x31*xp2**2 - 4*x11*x32*x41*xp2 + 2*x11*x32*x42*xp1 + 2*x11*x32*xp1*xp2 + 2*x11*x41*x42*xp2 + 2*x11*x41*xp2**2 - 2*x11*x42**2*xp1 - 2*x11*x42*xp1*xp2 + x12**2*x41**2 - 2*x12**2*x41*xp1 + x12**2*xp1**2 + 4*x12*x21*x31*x42 - 4*x12*x21*x31*xp2 - 2*x12*x21*x32*x41 + 2*x12*x21*x32*xp1 + 2*x12*x21*x41*xp2 - 4*x12*x21*x42*xp1 + 2*x12*x21*xp1*xp2 - 2*x12*x22*x31*x41 + 2*x12*x22*x31*xp1 + 2*x12*x22*x41*xp1 - 2*x12*x22*xp1**2 + 2*x12*x31*x41*xp2 - 4*x12*x31*x42*xp1 + 2*x12*x31*xp1*xp2 + 2*x12*x32*x41*xp1 - 2*x12*x32*xp1**2 - 2*x12*x41**2*xp2 + 2*x12*x41*x42*xp1 - 2*x12*x41*xp1*xp2 + 2*x12*x42*xp1**2 + x21**2*x32**2 - 2*x21**2*x32*xp2 + x21**2*xp2**2 - 2*x21*x22*x31*x32 + 2*x21*x22*x31*xp2 + 2*x21*x22*x32*xp1 - 2*x21*x22*xp1*xp2 + 2*x21*x31*x32*xp2 - 4*x21*x31*x42*xp2 + 2*x21*x31*xp2**2 - 2*x21*x32**2*xp1 + 2*x21*x32*x41*xp2 + 2*x21*x32*x42*xp1 - 2*x21*x32*xp1*xp2 - 2*x21*x41*xp2**2 + 2*x21*x42*xp1*xp2 + x22**2*x31**2 - 2*x22**2*x31*xp1 + x22**2*xp1**2 - 2*x22*x31**2*xp2 + 2*x22*x31*x32*xp1 + 2*x22*x31*x41*xp2 + 2*x22*x31*x42*xp1 - 2*x22*x31*xp1*xp2 - 4*x22*x32*x41*xp1 + 2*x22*x32*xp1**2 + 2*x22*x41*xp1*xp2 - 2*x22*x42*xp1**2 + x31**2*xp2**2 - 2*x31*x32*xp1*xp2 - 2*x31*x41*xp2**2 + 2*x31*x42*xp1*xp2 + x32**2*xp1**2 + 2*x32*x41*xp1*xp2 - 2*x32*x42*xp1**2 + x41**2*xp2**2 - 2*x41*x42*xp1*xp2 + x42**2*xp1**2))/(2*(x11*x32 - x11*x42 - x12*x31 + x12*x41 - x21*x32 + x21*x42 + x22*x31 - x22*x41))
  
  print "xi1=",expr
  xi1_of_xi3 = expr.subs([
    (x11,(1-xi3)*x11 + xi3*x51), (x21,(1-xi3)*x21 + xi3*x61), (x31,(1-xi3)*x31 + xi3*x71), (x41,(1-xi3)*x41 + xi3*x81),
    (x12,(1-xi3)*x12 + xi3*x52), (x22,(1-xi3)*x22 + xi3*x62), (x32,(1-xi3)*x32 + xi3*x72), (x42,(1-xi3)*x42 + xi3*x82)])
  
  print "-"
  xi1_of_xi3b = expr.subs([
    (x11,(1-xi3)*x11 + xi3*x51), (x21,(1-xi3)*x21 + xi3*x61), (x31,(1-xi3)*x31 + xi3*x71), (x41,(1-xi3)*x41 + xi3*x81),
    (x12,(1-xi3)*x13 + xi3*x53), (x22,(1-xi3)*x23 + xi3*x63), (x32,(1-xi3)*x33 + xi3*x73), (x42,(1-xi3)*x43 + xi3*x83)])
  print "-"
  
  # xi2
  expr = (2*x11*x22 - x11*x42 - 2*x12*x21 + x12*x41 + x21*x32 - x22*x31 + xp1*(x12 - x22 - x32 + x42) + xp2*(-x11 + x21 + x31 - x41) - sqrt(x11**2*x42**2 - 2*x11**2*x42*xp2 + x11**2*xp2**2 - 2*x11*x12*x41*x42 + 2*x11*x12*x41*xp2 + 2*x11*x12*x42*xp1 - 2*x11*x12*xp1*xp2 - 2*x11*x21*x32*x42 + 2*x11*x21*x32*xp2 + 2*x11*x21*x42*xp2 - 2*x11*x21*xp2**2 - 2*x11*x22*x31*x42 + 2*x11*x22*x31*xp2 + 4*x11*x22*x32*x41 - 4*x11*x22*x32*xp1 - 4*x11*x22*x41*xp2 + 2*x11*x22*x42*xp1 + 2*x11*x22*xp1*xp2 + 2*x11*x31*x42*xp2 - 2*x11*x31*xp2**2 - 4*x11*x32*x41*xp2 + 2*x11*x32*x42*xp1 + 2*x11*x32*xp1*xp2 + 2*x11*x41*x42*xp2 + 2*x11*x41*xp2**2 - 2*x11*x42**2*xp1 - 2*x11*x42*xp1*xp2 + x12**2*x41**2 - 2*x12**2*x41*xp1 + x12**2*xp1**2 + 4*x12*x21*x31*x42 - 4*x12*x21*x31*xp2 - 2*x12*x21*x32*x41 + 2*x12*x21*x32*xp1 + 2*x12*x21*x41*xp2 - 4*x12*x21*x42*xp1 + 2*x12*x21*xp1*xp2 - 2*x12*x22*x31*x41 + 2*x12*x22*x31*xp1 + 2*x12*x22*x41*xp1 - 2*x12*x22*xp1**2 + 2*x12*x31*x41*xp2 - 4*x12*x31*x42*xp1 + 2*x12*x31*xp1*xp2 + 2*x12*x32*x41*xp1 - 2*x12*x32*xp1**2 - 2*x12*x41**2*xp2 + 2*x12*x41*x42*xp1 - 2*x12*x41*xp1*xp2 + 2*x12*x42*xp1**2 + x21**2*x32**2 - 2*x21**2*x32*xp2 + x21**2*xp2**2 - 2*x21*x22*x31*x32 + 2*x21*x22*x31*xp2 + 2*x21*x22*x32*xp1 - 2*x21*x22*xp1*xp2 + 2*x21*x31*x32*xp2 - 4*x21*x31*x42*xp2 + 2*x21*x31*xp2**2 - 2*x21*x32**2*xp1 + 2*x21*x32*x41*xp2 + 2*x21*x32*x42*xp1 - 2*x21*x32*xp1*xp2 - 2*x21*x41*xp2**2 + 2*x21*x42*xp1*xp2 + x22**2*x31**2 - 2*x22**2*x31*xp1 + x22**2*xp1**2 - 2*x22*x31**2*xp2 + 2*x22*x31*x32*xp1 + 2*x22*x31*x41*xp2 + 2*x22*x31*x42*xp1 - 2*x22*x31*xp1*xp2 - 4*x22*x32*x41*xp1 + 2*x22*x32*xp1**2 + 2*x22*x41*xp1*xp2 - 2*x22*x42*xp1**2 + x31**2*xp2**2 - 2*x31*x32*xp1*xp2 - 2*x31*x41*xp2**2 + 2*x31*x42*xp1*xp2 + x32**2*xp1**2 + 2*x32*x41*xp1*xp2 - 2*x32*x42*xp1**2 + x41**2*xp2**2 - 2*x41*x42*xp1*xp2 + x42**2*xp1**2))/(2*(x11*x22 - x11*x42 - x12*x21 + x12*x41 + x21*x32 - x22*x31 + x31*x42 - x32*x41))
  xi1_of_xi2 = expr.subs([
    (x11,(1-xi3)*x11 + xi3*x51), (x21,(1-xi3)*x21 + xi3*x61), (x31,(1-xi3)*x31 + xi3*x71), (x41,(1-xi3)*x41 + xi3*x81),
    (x12,(1-xi3)*x12 + xi3*x52), (x22,(1-xi3)*x22 + xi3*x62), (x32,(1-xi3)*x32 + xi3*x72), (x42,(1-xi3)*x42 + xi3*x82)])
  
  print "-"
  xi1_of_xi3 = collect(xi1_of_xi3,xi3)
  print "-"
  xi1_of_xi3b = collect(xi1_of_xi3b,xi3)
  print "-"
  xi1_of_xi2 = collect(xi1_of_xi2,xi2)
  
  print "xi1_of_xi3=",xi1_of_xi3
  print ""
  print "xi1_of_xi3b=",xi1_of_xi3b
  print ""
  print "xi1_of_xi2=",xi1_of_xi2
  print ""
  
  print "solve for xi3:"
  xi3sol = solve(xi1_of_xi3,xi3)
  print "xi3=",xi3sol
  print ""
  
  print "substitute xi3 in xi1b"
  xi1sol = xi1_of_xi3b.subs(xi3,xi3sol)
  print "xi1=",xi1sol
  print ""
  
  print "solve for xi2:"
  xi2sol = solve(xi1_of_xi2,xi2)
  print "xi2=",xi2sol
  print ""
  
if False:
  Phi1 = (1-xi1)*(1-xi2)*(1-xi3)*x11 + xi1*(1-xi2)*(1-xi3)*x21 + (1-xi1)*xi2*(1-xi3)*x31 + xi1*xi2*(1-xi3)*x41 + (1-xi1)*(1-xi2)*xi3*x51 + xi1*(1-xi2)*xi3*x61 + (1-xi1)*xi2*xi3*x71 + xi1*xi2*xi3*x81
  Phi2 = (1-xi1)*(1-xi2)*(1-xi3)*x12 + xi1*(1-xi2)*(1-xi3)*x22 + (1-xi1)*xi2*(1-xi3)*x32 + xi1*xi2*(1-xi3)*x42 + (1-xi1)*(1-xi2)*xi3*x52 + xi1*(1-xi2)*xi3*x62 + (1-xi1)*xi2*xi3*x72 + xi1*xi2*xi3*x82
  Phi3 = (1-xi1)*(1-xi2)*(1-xi3)*x13 + xi1*(1-xi2)*(1-xi3)*x23 + (1-xi1)*xi2*(1-xi3)*x33 + xi1*xi2*(1-xi3)*x43 + (1-xi1)*(1-xi2)*xi3*x53 + xi1*(1-xi2)*xi3*x63 + (1-xi1)*xi2*xi3*x73 + xi1*xi2*xi3*x83

  #Phi1 = xi1+xi2*xi3
  #Phi2 = xi1 + x11 + xi2+x12 + xi3
  #Phi3 = x13+(xi2+xi3) + 2*xi1
  
  solution = solve([Phi1 - xp1, Phi2 - xp2, Phi3 - xp3], [xi1,xi2,xi3])
  print solution
  print ""
  
  sys.exit(0)
  
  solution13 = solve([Phi1 - xp1, Phi2 - xp2], [xi1,xi3])
  print solution
  print ""
  
  sys.exit(0)
  print "3D linear"
  #solution = solve([Phi1 - xp1, Phi2 - xp2, Phi3 - xp3], [xi1,xi2,xi3])
  #print len(solution)," solution(s)"
  #print solution
  #print "----------"
  #solution = simplify(solution)
  #print solution
  #print ""
  
  xi1sol = solve([xp1-Phi1], [xi1])[xi1]
  print "xi1: ", xi1sol
  # xi1sol(xi2,xi3)
  
  print "sustitute Phi2sol = Phi2.subs(xi1,xi1sol)"
  Phi2sol = (t1 + xi2*t2 + xi2**2*t3) / (t4 + xi2*t5)
  # Phi2(xi1,xi2,xi3
  # Phi2sol(xi2,xi3)
  
  t1 = -x11*x12*xi3*xia1 + x11*x12*xia1 + x11*x22*xi3**2 + x11*x22*xi3*xia1 - 2*x11*x22*xi3 - x11*x22*xia1 + x11*x22 + x11*x52*xi3*xia1 - x11*x62*xi3**2 - x11*x62*xi3*xia1 + x11*x62*xi3 - x12*x21*xi3**2 + x12*x21*xi3*xia1 + x12*x21*xi3*xia4 + x12*x21*xi3*xia7 + x12*x21*xi3 - x12*x21*xia1 - x12*x21*xia4 - x12*x21*xia7 - x12*x41*xi3*xia4 + x12*x41*xia4 + x12*x61*xi3**2 - x12*x61*xi3 - x12*x71*xi3*xia7 + x12*x71*xia7 - x21*x22*xi3*xia1 - x21*x22*xi3*xia4 - x21*x22*xi3*xia7 + x21*x22*xi3 + x21*x22*xia1 + x21*x22*xia4 + x21*x22*xia7 - x21*x22 + x21*x52*xi3**2 - x21*x52*xi3*xia1 - x21*x52*xi3*xia4 - x21*x52*xi3*xia7 + x21*x62*xi3*xia1 + x21*x62*xi3*xia4 + x21*x62*xi3*xia7 - x21*x62*xi3 + x22*x41*xi3*xia4 - x22*x41*xia4 - x22*x51*xi3**2 + x22*x51*xi3 + x22*x71*xi3*xia7 - x22*x71*xia7 + x41*x52*xi3*xia4 - x41*x62*xi3*xia4 + x51*x62*xi3**2 - x52*x61*xi3**2 + x52*x71*xi3*xia7 - x62*x71*xi3*xia7 
  t2 = x11*x12*xi3*xia1 - x11*x12*xia1 - 2*x11*x22*xi3**2 - x11*x22*xi3*xia1 + 4*x11*x22*xi3 + x11*x22*xia1 - 2*x11*x22 - x11*x32*xi3*xia1 + x11*x32*xia1 + x11*x42*xi3**2 + x11*x42*xi3*xia1 - 2*x11*x42*xi3 - x11*x42*xia1 + x11*x42 - x11*x52*xi3*xia1 + 2*x11*x62*xi3**2 + x11*x62*xi3*xia1 - 2*x11*x62*xi3 + x11*x72*xi3*xia1 - x11*x82*xi3**2 - x11*x82*xi3*xia1 + x11*x82*xi3 + 2*x12*x21*xi3**2 - x12*x21*xi3*xia1 - x12*x21*xi3*xia4 - x12*x21*xi3*xia7 - 3*x12*x21*xi3 + x12*x21*xia1 + x12*x21*xia4 + x12*x21*xia7 + x12*x21 - x12*x41*xi3**2 + x12*x41*xi3*xia4 + 2*x12*x41*xi3 - x12*x41*xia4 - x12*x41 - 2*x12*x61*xi3**2 + 2*x12*x61*xi3 + x12*x71*xi3*xia7 - x12*x71*xia7 + x12*x81*xi3**2 - x12*x81*xi3 + x21*x22*xi3*xia1 + x21*x22*xi3*xia4 + x21*x22*xi3*xia7 - x21*x22*xi3 - x21*x22*xia1 - x21*x22*xia4 - x21*x22*xia7 + x21*x22 - x21*x32*xi3**2 + x21*x32*xi3*xia1 + x21*x32*xi3*xia4 + x21*x32*xi3*xia7 + x21*x32*xi3 - x21*x32*xia1 - x21*x32*xia4 - x21*x32*xia7 - x21*x42*xi3*xia1 - x21*x42*xi3*xia4 - x21*x42*xi3*xia7 + x21*x42*xi3 + x21*x42*xia1 + x21*x42*xia4 + x21*x42*xia7 - x21*x42 - 2*x21*x52*xi3**2 + x21*x52*xi3*xia1 + x21*x52*xi3*xia4 + x21*x52*xi3*xia7 + x21*x52*xi3 - x21*x62*xi3*xia1 - x21*x62*xi3*xia4 - x21*x62*xi3*xia7 + x21*x62*xi3 + x21*x72*xi3**2 - x21*x72*xi3*xia1 - x21*x72*xi3*xia4 - x21*x72*xi3*xia7 + x21*x82*xi3*xia1 + x21*x82*xi3*xia4 + x21*x82*xi3*xia7 - x21*x82*xi3 + x22*x31*xi3**2 - 2*x22*x31*xi3 + x22*x31 - x22*x41*xi3*xia4 + x22*x41*xia4 + 2*x22*x51*xi3**2 - 2*x22*x51*xi3 - x22*x71*xi3**2 - x22*x71*xi3*xia7 + x22*x71*xi3 + x22*x71*xia7 - x31*x62*xi3**2 + x31*x62*xi3 - x32*x41*xi3*xia4 + x32*x41*xia4 + x32*x61*xi3**2 - x32*x61*xi3 - x32*x71*xi3*xia7 + x32*x71*xia7 + x41*x42*xi3*xia4 - x41*x42*xia4 + x41*x52*xi3**2 - x41*x52*xi3*xia4 - x41*x52*xi3 + x41*x62*xi3*xia4 + x41*x72*xi3*xia4 - x41*x82*xi3*xia4 - x42*x51*xi3**2 + x42*x51*xi3 + x42*x71*xi3*xia7 - x42*x71*xia7 - 2*x51*x62*xi3**2 + x51*x82*xi3**2 + 2*x52*x61*xi3**2 - x52*x71*xi3*xia7 - x52*x81*xi3**2 - x61*x72*xi3**2 + x62*x71*xi3**2 + x62*x71*xi3*xia7 + x71*x72*xi3*xia7 - x71*x82*xi3*xia7
  t3 = x11*x22*xi3**2 - 2*x11*x22*xi3 + x11*x22 - x11*x42*xi3**2 + 2*x11*x42*xi3 - x11*x42 - x11*x62*xi3**2 + x11*x62*xi3 + x11*x82*xi3**2 - x11*x82*xi3 - x12*x21*xi3**2 + 2*x12*x21*xi3 - x12*x21 + x12*x41*xi3**2 - 2*x12*x41*xi3 + x12*x41 + x12*x61*xi3**2 - x12*x61*xi3 - x12*x81*xi3**2 + x12*x81*xi3 + x21*x32*xi3**2 - 2*x21*x32*xi3 + x21*x32 + x21*x52*xi3**2 - x21*x52*xi3 - x21*x72*xi3**2 + x21*x72*xi3 - x22*x31*xi3**2 + 2*x22*x31*xi3 - x22*x31 - x22*x51*xi3**2 + x22*x51*xi3 + x22*x71*xi3**2 - x22*x71*xi3 + x31*x42*xi3**2 - 2*x31*x42*xi3 + x31*x42 + x31*x62*xi3**2 - x31*x62*xi3 - x31*x82*xi3**2 + x31*x82*xi3 - x32*x41*xi3**2 + 2*x32*x41*xi3 - x32*x41 - x32*x61*xi3**2 + x32*x61*xi3 + x32*x81*xi3**2 - x32*x81*xi3 - x41*x52*xi3**2 + x41*x52*xi3 + x41*x72*xi3**2 - x41*x72*xi3 + x42*x51*xi3**2 - x42*x51*xi3 - x42*x71*xi3**2 + x42*x71*xi3 + x51*x62*xi3**2 - x51*x82*xi3**2 - x52*x61*xi3**2 + x52*x81*xi3**2 + x61*x72*xi3**2 - x62*x71*xi3**2 + x71*x82*xi3**2 - x72*x81*xi3**2
  t4 = -x11*xi3 + x11 + x21*xi3 - x21 + x51*xi3 - x61*xi3
  t5 = x11*xi3 - x11 - x21*xi3 + x21 - x31*xi3 + x31 + x41*xi3 - x41 - x51*xi3 + x61*xi3 + x71*xi3 - x81*xi3
  
  xi2sol = solve([xp2-Phi2sol], [xi2])[0]
  
  xi2_a = (-t2 + t5*xp2 - sqrt(-4*t1*t3 + t2**2 - 2*t2*t5*xp2 + 4*t3*t4*xp2 + t5**2*xp2**2))/(2*t3)
  xi2_b = (-t2 + t5*xp2 + sqrt(-4*t1*t3 + t2**2 - 2*t2*t5*xp2 + 4*t3*t4*xp2 + t5**2*xp2**2))/(2*t3)
  
  # xi2sol(xi3)
  
  
  

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

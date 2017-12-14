#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
from sympy.matrices import *
from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False, no_global=True)

m11,m12,m13 = symbols('m11,m12,m13')
m21,m22,m23 = symbols('m21,m22,m23')
m31,m32,m33 = symbols('m31,m32,m33')

# compute T
M=Matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])
print M
print "inv:",simplify(M.inv())
# inv = 1/det * adj

print "adj:",simplify(M.adjugate())
 
print "det:",simplify(M.det())
print "T: 1/(det**2)*",simplify(M.adjugate()*M.adjugate().transpose())
 
# 3D
# compute T
l1,l2 = symbols('l1,l2')
M2 = Matrix([[l1, -l1*cos(a)/sin(a)],[0,l2/sin(a)]])
M2inv = M2.inv()
print "M2:",M2
print "inv:",M2inv
print "invT:",M2inv.transpose()
print "T:",M2inv*M2inv.transpose()


# compute v1*S*v2, S is a symmetric matrix
v11,v12,v13 = symbols('v11,v12,v13')
v21,v22,v23 = symbols('v21,v22,v23')
S=Matrix([[m11,m12,m13],[m12,m22,m23],[m13,m23,m33]])
v1=Matrix([[v11],[v12],[v13]])
v2=Matrix([[v21],[v22],[v23]])
print simplify(v1.transpose()*S*v2)
# v21*(m11*v11 + m12*v12 + m13*v13) + v22*(m12*v11 + m22*v12 + m23*v13) + v23*(m13*v11 + m23*v12 + m33*v13)

# 2D
M = Matrix([[m11,m12],[m21,m22],[m31,m32]])
print "J2:",sqrt(det(M.transpose()*M))

S=Matrix([[m11,m12],[m12,m22]])
v1=Matrix([[v11],[v12]])
v2=Matrix([[v21],[v22]])
print simplify(v1.transpose()*S*v2)

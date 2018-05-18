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

a11,a12,a13 = symbols('a11,a12,a13')
a21,a22,a23 = symbols('a21,a22,a23')
a31,a32,a33 = symbols('a31,a32,a33')

v11,v12,v13 = symbols('v11,v12,v13')
v21,v22,v23 = symbols('v21,v22,v23')

##################
# laplace operator

# compute v1^T * M * v2, M is a symmetric matrix, the integral transformation matrix
# 3D
M = Matrix([[m11,m12,m13],[m12,m22,m23],[m13,m23,m33]])
v1 = Matrix([[v11],[v12],[v13]])
v2 = Matrix([[v21],[v22],[v23]])
print ""
print simplify(v1.transpose()*M*v2)
# v21*(m11*v11 + m12*v12 + m13*v13) + v22*(m12*v11 + m22*v12 + m23*v13) + v23*(m13*v11 + m23*v12 + m33*v13)

# 2D
M = Matrix([[m11,m12],[m21,m22],[m31,m32]])
print ""
print "J2:",sqrt(det(M.transpose()*M))

M = Matrix([[m11,m12],[m12,m22]])
v1 = Matrix([[v11],[v12]])
v2 = Matrix([[v21],[v22]])
print ""
print simplify(v1.transpose()*M*v2)

##############################
# generalized laplace operator

# compute A * v1^T * M * v2, M is a symmetric matrix, the integral transformation matrix, A is the conductivity tensor
# 3D
M = Matrix([[m11,m12,m13],[m12,m22,m23],[m13,m23,m33]])
A = Matrix([[a11,a12,a13],[a12,a22,a23],[a13,a23,a33]])
v1 = Matrix([[v11],[v12],[v13]])
v2 = Matrix([[v21],[v22],[v23]])
Av1 = A*v1
print ""
print simplify(Av1.transpose()*M*v2)

# 2D
M = Matrix([[m11,m12],[m21,m22],[m31,m32]])

M = Matrix([[m11,m12],[m12,m22]])
A = Matrix([[a11,a12],[a12,a22]])
v1 = Matrix([[v11],[v12]])
v2 = Matrix([[v21],[v22]])
Av1 = A*v1
print ""
print simplify(Av1.transpose()*M*v2)

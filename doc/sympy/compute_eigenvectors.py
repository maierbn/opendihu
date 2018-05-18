#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *

# symmetric matrix A
a11,a12,a13,a21,a22,a23,a31,a32,a33 = symbols('m11,m12,m13,m21,m22,m23,m31,m32,m33')

# inverse
matrix = Matrix([
  [a11, a21, a31],
  [a21, a22, a32],
  [a31, a32, a33]
])
det = det(matrix)
adj = matrix.inv()*det
inverse = matrix.inv()
print "det(matrix):",det
print ""
print "adj(matrix):",simplify(adj)
print ""
print "inv(matrix):",inverse

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

# compute det(M) * M^{-T} = det(M) * inv^T = det(M) * 1/det * cof^TT = cof
# compute inverse of M
M=Matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])
print "matrix:"
print M
print ""
print "inv:",simplify(M.inv())
# inv = 1/det * adj = 1/det * cof^T

print "det:",simplify(M.det())
print ""
print "adj:",simplify(M.adjugate())
print ""
print "cof:",simplify(M.cofactorMatrix())
 


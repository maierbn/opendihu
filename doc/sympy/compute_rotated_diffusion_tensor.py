#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
from symby.vector import CoordSys3D
N = CoordSys3D('N')

# direction vector
# rotation from a to b
a1,a2,a3 = symbols('a1,a2,a3')
b1,b2,b3 = symbols('b1,b2,b3')

a = a1*N.i + a2*N.j + a3*N.k
b = b1*N.i + b2*N.j + b3*N.k

u = (a.dot(b))*a / norm((a.dot(b))*a).norm()
v = (b - a.dot(b))*a) / (b - (a.dot(b))*a).norm()
w = b.cross(a)

f = Matrix([
  [u[0], v[0], w[0]],
  [u[1], v[1], w[1]],
  [u[2], v[2], w[2]]
]).inv()

g = Matrix([
  [a.dot(b),          -(a.cross(b)).norm()), 0],
  [a.cross(b).norm(), -a.dot(b),             0],
  [0,                 0,                     1]
])

uu = f.inv()*g*f

print ""
print "f:"
print f
print ""
print "g:"
print g
print ""
print "U:"
print uu
print ""
print "Ua:"
print u.dot(a)

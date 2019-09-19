#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
from sympy.vector import CoordSys3D
N = CoordSys3D('N')

# generate rotation matrix uu, such that  uu*a = b, a = [1,0,0]
# rotation is done by an axis w = b x a, in the a-b-plane
# the rotate coordinate system is [u,v,w] with u = (a•b)*a, v = b - (a•b)*a, w = b x a (normalized each)
# the rotation matrix g is [[cos(phi), -sin(phi), 0], [sin(phi), cos(phi), 0], [0, 0, 1]], where phi = /_(a,b)


# direction vector
# rotation from a to b
a1,a2,a3 = symbols('a1,a2,a3', positive=True)
b1,b2,b3 = symbols('b1,b2,b3', positive=True)

#a = Matrix([a1,a2,a3])
a = Matrix([1,0,0])

b = Matrix([b1,b2,b3])

uhelper = (a.dot(b))*a
u = uhelper / uhelper.norm()
print ""
print "u=",u


vhelper = b - a.dot(b)*a
v = vhelper / vhelper.norm()
print ""
print "v=",v

w = b.cross(a)
print ""
print "w=",w

finv = Matrix([
  [u[0], v[0], w[0]],
  [u[1], v[1], w[1]],
  [u[2], v[2], w[2]]
])
f = factor(finv.inv())

print ""
print "f=",f

g = factor(Matrix([
  [a.dot(b),          -(a.cross(b)).norm(),  0],
  [a.cross(b).norm(), -a.dot(b),             0],
  [0,                 0,                     1]
]))

print ""
print "g=",g

uu = factor(f.inv()*g*f)

print ""
print "f:"
print f
print ""
print "g:"
print g

print ""
print "transformation matrix U 3D:"
print uu
print ""
print "inverse:"
print simplify(factor(uu.inv()))
print "det=",simplify(factor(uu.det()))

print ""
print "transformation matrix U 2D:"
print uu.subs(b3,0)
print ""
print "inverse"
print factor(uu.subs(b3,0).inv())

print ""
print "Ua:"
print uu.dot(a)
print ""
print factor(uu.dot(a))

# diffusion tensor
import numpy as np
uu_num = uu.subs([(b1,0.), (b2,0.), (b3,5.)])

print "b=",b.subs([(b1,0.), (b2,0.), (b3,5.)])

base_matrix = Matrix([[2,0,0], [0,1,0], [0,0,1]])
print("rotation matrix: ",uu_num)
print("rotation_matrix inverse: ",uu_num.inv())

sigma = uu_num.inv()*base_matrix*uu_num
print "rotInv*matrix:",uu_num.inv()*base_matrix
print "sigma=",sigma
eigenvectors = sigma.eigenvects(simplify=True)
for eigenvector in eigenvectors:
  print ""
  print "eigenval: ",eigenvector[0]
  print "multiplicity: ",eigenvector[1]
  print "eigenvector: ",eigenvector[2]


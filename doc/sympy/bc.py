#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# weak form Dirichlet BC, Δu = 0, 5 elements, length = 4, linear ansatz functions, u(0) = 1, u(5) = 0, solution: [1, 0.8, 0.6, 0.4, 0.2, 0.0]
#

import numpy as np

m = np.array([[-5./4, 5./4, 0, 0, 0, 0], [5./4, -5./2, 5./4, 0, 0, 0], [0, 5./4, -5./2, 5./4, 0, 0], [0, 0, 5./4, -5./2, 5./4, 0], [0, 0, 0, 5./4, -5./2, 5./4], [0, 0, 0, 0, 5./4, -5./4]])
rhs = np.array([[-5./4],[5./4],[0.],[0],[0],[0]])
rhs = np.array([[0],[0],[0],[0],[0],[0]])

rhs2 = np.zeros([6,1])
m2 = np.array(m)

u0 = 1.0
u5 = 5.0

print "u0:",u0
print "u5:",u5
rhs2[0,0] = u0
rhs2[5,0] = u5

n = 6
for i in range(1,n-1):
  m2[0,i] = 0
  m2[5,i] = 0
  m2[i,0] = 0
  m2[i,5] = 0
  
  rhs2[i,0] -= m[i,0]*u0
  rhs2[i,0] -= m[i,5]*u5
  print "m[{},0]={}, set rhs2[{},0] to {}={}".format(i,m[i,0],i,-m[i,0]*u0,rhs2[i,0])
  
m2[0,0] = 1.0
m2[5,5] = 1.0
  
mm = np.zeros([6,2])
mm = m[:,[0,5]]
print "mm=",mm
rhs3 = -mm.dot(np.array([[u0],[u5]]))
rhs3[0] = u0
rhs3[5] = u5

print "rhs3=",rhs3
print rhs
print m
print ""
print rhs2
print m2
u2 = np.linalg.solve(m2,rhs2)
print u2
print ""
print ""
print "probe m*u2"
print m.dot(u2)

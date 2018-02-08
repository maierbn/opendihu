#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
import numpy as np
import sys

## 2D case: analytically compute sqrt(Gram(x1,x2))
x11,x12,x13,x21,x22,x23,x31,x32,x33,x41,x42,x43 = symbols('x11,x12,x13,x21,x22,x23,x31,x32,x33,x41,x42,x43')
xi1,xi2 = symbols('xi1, xi2')

def Phi(xi1,xi2):
  p1 = (1-xi1)*(1-xi2)*x11 + xi1*(1-xi2)*x21 + (1-xi1)*xi2*x31 + xi1*xi2*x41
  p2 = (1-xi1)*(1-xi2)*x12 + xi1*(1-xi2)*x22 + (1-xi1)*xi2*x32 + xi1*xi2*x42
  p3 = (1-xi1)*(1-xi2)*x13 + xi1*(1-xi2)*x23 + (1-xi1)*xi2*x33 + xi1*xi2*x43
  
  return [p1,p2,p3]

dphi_dxi1 = np.array([diff(Phi(xi1,xi2)[0],xi1), diff(Phi(xi1,xi2)[1],xi1), diff(Phi(xi1,xi2)[2],xi1)])
dphi_dxi2 = np.array([diff(Phi(xi1,xi2)[0],xi2), diff(Phi(xi1,xi2)[1],xi2), diff(Phi(xi1,xi2)[2],xi2)])

term = np.dot(dphi_dxi1,dphi_dxi1)*np.dot(dphi_dxi2,dphi_dxi2) - np.dot(dphi_dxi1,dphi_dxi2)**2
print term
j = simplify(sqrt(term))

print "2D J(xi):"
print j

## 3D case: analytically compute |det(DPhi(xi))| 
x11,x12,x13,x21,x22,x23,x31,x32,x33,x41,x42,x43 = symbols('x11,x12,x13,x21,x22,x23,x31,x32,x33,x41,x42,x43')
xi1,xi2,xi3 = symbols('xi1,xi2,xi3')

def Phi(xi1,xi2,xi3):
  p1 = (1-xi1)*(1-xi2)*(1-xi3)*x11 + xi1*(1-xi2)*(1-xi3)*x21 + (1-xi1)*xi2*(1-xi3)*x31 + xi1*xi2*(1-xi3)*x41 + (1-xi1)*(1-xi2)*xi3*x11 + xi1*(1-xi2)*xi3*x21 + (1-xi1)*xi2*xi3)*x31 + xi1*xi2*xi3)*x41
  p2 = (1-xi1)*(1-xi2)*(1-xi3)*x12 + xi1*(1-xi2)*(1-xi3)*x22 + (1-xi1)*xi2*(1-xi3)*x32 + xi1*xi2*(1-xi3)*x42 + (1-xi1)*(1-xi2)*xi3*x12 + xi1*(1-xi2)*xi3*x22 + (1-xi1)*xi2*xi3)*x32 + xi1*xi2*xi3)*x42
  p3 = (1-xi1)*(1-xi2)*(1-xi3)*x13 + xi1*(1-xi2)*(1-xi3)*x23 + (1-xi1)*xi2*(1-xi3)*x33 + xi1*xi2*(1-xi3)*x43 + (1-xi1)*(1-xi2)*xi3*x13 + xi1*(1-xi2)*xi3*x23 + (1-xi1)*xi2*xi3)*x33 + xi1*xi2*xi3)*x43
  
  return [p1,p2,p3]

matrix = Matrix([
  [diff(Phi(xi1,xi2,xi3)[0],xi1), diff(Phi(xi1,xi2,xi3)[0],xi2), diff(Phi(xi1,xi2,xi3)[0],xi3)],
  [diff(Phi(xi1,xi2,xi3)[1],xi1), diff(Phi(xi1,xi2,xi3)[1],xi2), diff(Phi(xi1,xi2,xi3)[1],xi3)]
  [diff(Phi(xi1,xi2,xi3)[2],xi1), diff(Phi(xi1,xi2,xi3)[2],xi2), diff(Phi(xi1,xi2,xi3)[2],xi3)]
])
term = det(matrix)
print "det(DPhi(xi)):"
print collect(term,(xi1,xi2))

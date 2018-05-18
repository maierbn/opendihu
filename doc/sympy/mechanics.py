#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *
D = 2
C11, C21, C12, C22 = symbols('C11, C21, C12, C22')
c0, c1 = symbols('c0, c1')
p,l = symbols('p,l')


C = Matrix([[C11, C21],[C12, C22]])

# set C for problem
C11 = l**2
C22 = l**(-2)
C12 = 0
C21 = 0

C = Matrix([[C11, C21],[C12, C22]])

print "C: ", C
print ""


# derive S analytically
def deriveS(C):
  global Svol, Siso
  
  I = Matrix([[1,0],[0,1]])
  Cinv = C.inv()

  Cbar = C
  gammabar1 = 2*c0
  gammabar2 = 0
  Sbar = gammabar1*I + gammabar2*Cbar


  CSbar = C[0] * Sbar[0]
  for i in range(1,D*D):
    CSbar += C[i] * Sbar[i]
  print "CSbar: ", CSbar
  print ""

  Siso = Sbar - 1./2*Cinv*CSbar
  Siso = simplify(Siso)
  Svol = p * Cinv
  S = Siso + Svol

  print "Siso: ", Siso
  print ""
  print "Svol: ", Svol
  print ""
  
  return S
  
S = deriveS(C)
print "S: ", S
print ""

# get analytic value for p
p_analytic = solve(S[D*D-1], p, simplify=True, positive=True)[0]
print "p analytic: ", p_analytic
print ""

# substitute p in Svol
Svol_analytic = simplify(Svol.subs([(p, p_analytic)]))
print "Svol analytic", Svol_analytic
print ""



# substitute p in S
S_analytic = simplify(S.subs([(p, p_analytic)]))
print "S analytic", S_analytic
print ""


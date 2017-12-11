#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.abc import *

dPhi1_dxi1, dPhi2_dxi1, dPhi3_dxi1 = symbols('dPhi1_dxi1, dPhi2_dxi1, dPhi3_dxi1')
dPhi1_dxi2, dPhi2_dxi2, dPhi3_dxi2 = symbols('dPhi1_dxi2, dPhi2_dxi2, dPhi3_dxi2')

zeta1 = Matrix([dPhi1_dxi1, dPhi2_dxi1, dPhi3_dxi1])
zetah = Matrix([dPhi1_dxi2, dPhi2_dxi2, dPhi3_dxi2])
zeta3 = zeta1.cross(zetah)
zeta2 = zeta3.cross(zeta1)

print "zeta1: ",zeta1
print "zetah: ",zetah
print "zeta3: ",zeta3
print "zeta2: ",simplify(zeta2)

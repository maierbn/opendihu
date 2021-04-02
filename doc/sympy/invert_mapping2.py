#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys

def get_xi_2d(quadrilateral,p):
  [p0, p1, p2, p3] = quadrilateral

  xp1 = p[0]
  xp2 = p[1]
  xp3 = p[2]
  
  x11 = p0[0]
  x12 = p0[1]
  x13 = p0[2]
  x21 = p1[0]
  x22 = p1[1]
  x23 = p1[2]
  x31 = p2[0]
  x32 = p2[1]
  x33 = p2[2]
  x41 = p3[0]
  x42 = p3[1]
  x43 = p3[2]
  
  xi1 = (2*x11*x32 - x11*x42 - 2*x12*x31 + x12*x41 - x21*x32 + x22*x31 + xp1*(x12 - x22 - x32 + x42) + xp2*(-x11 + x21 + x31 - x41) + np.sqrt(x11**2*x42**2 - 2*x11**2*x42*xp2 + x11**2*xp2**2 - 2*x11*x12*x41*x42 + 2*x11*x12*x41*xp2 + 2*x11*x12*x42*xp1 - 2*x11*x12*xp1*xp2 - 2*x11*x21*x32*x42 + 2*x11*x21*x32*xp2 + 2*x11*x21*x42*xp2 - 2*x11*x21*xp2**2 - 2*x11*x22*x31*x42 + 2*x11*x22*x31*xp2 + 4*x11*x22*x32*x41 - 4*x11*x22*x32*xp1 - 4*x11*x22*x41*xp2 + 2*x11*x22*x42*xp1 + 2*x11*x22*xp1*xp2 + 2*x11*x31*x42*xp2 - 2*x11*x31*xp2**2 - 4*x11*x32*x41*xp2 + 2*x11*x32*x42*xp1 + 2*x11*x32*xp1*xp2 + 2*x11*x41*x42*xp2 + 2*x11*x41*xp2**2 - 2*x11*x42**2*xp1 - 2*x11*x42*xp1*xp2 + x12**2*x41**2 - 2*x12**2*x41*xp1 + x12**2*xp1**2 + 4*x12*x21*x31*x42 - 4*x12*x21*x31*xp2 - 2*x12*x21*x32*x41 + 2*x12*x21*x32*xp1 + 2*x12*x21*x41*xp2 - 4*x12*x21*x42*xp1 + 2*x12*x21*xp1*xp2 - 2*x12*x22*x31*x41 + 2*x12*x22*x31*xp1 + 2*x12*x22*x41*xp1 - 2*x12*x22*xp1**2 + 2*x12*x31*x41*xp2 - 4*x12*x31*x42*xp1 + 2*x12*x31*xp1*xp2 + 2*x12*x32*x41*xp1 - 2*x12*x32*xp1**2 - 2*x12*x41**2*xp2 + 2*x12*x41*x42*xp1 - 2*x12*x41*xp1*xp2 + 2*x12*x42*xp1**2 + x21**2*x32**2 - 2*x21**2*x32*xp2 + x21**2*xp2**2 - 2*x21*x22*x31*x32 + 2*x21*x22*x31*xp2 + 2*x21*x22*x32*xp1 - 2*x21*x22*xp1*xp2 + 2*x21*x31*x32*xp2 - 4*x21*x31*x42*xp2 + 2*x21*x31*xp2**2 - 2*x21*x32**2*xp1 + 2*x21*x32*x41*xp2 + 2*x21*x32*x42*xp1 - 2*x21*x32*xp1*xp2 - 2*x21*x41*xp2**2 + 2*x21*x42*xp1*xp2 + x22**2*x31**2 - 2*x22**2*x31*xp1 + x22**2*xp1**2 - 2*x22*x31**2*xp2 + 2*x22*x31*x32*xp1 + 2*x22*x31*x41*xp2 + 2*x22*x31*x42*xp1 - 2*x22*x31*xp1*xp2 - 4*x22*x32*x41*xp1 + 2*x22*x32*xp1**2 + 2*x22*x41*xp1*xp2 - 2*x22*x42*xp1**2 + x31**2*xp2**2 - 2*x31*x32*xp1*xp2 - 2*x31*x41*xp2**2 + 2*x31*x42*xp1*xp2 + x32**2*xp1**2 + 2*x32*x41*xp1*xp2 - 2*x32*x42*xp1**2 + x41**2*xp2**2 - 2*x41*x42*xp1*xp2 + x42**2*xp1**2))/(2*(x11*x32 - x11*x42 - x12*x31 + x12*x41 - x21*x32 + x21*x42 + x22*x31 - x22*x41))
  
  xi2a = (2*x11*x22 - x11*x42 - 2*x12*x21 + x12*x41 + x21*x32 - x22*x31 + xp1*(x12 - x22 - x32 + x42) + xp2*(-x11 + x21 + x31 - x41) - np.sqrt(x11**2*x42**2 - 2*x11**2*x42*xp2 + x11**2*xp2**2 - 2*x11*x12*x41*x42 + 2*x11*x12*x41*xp2 + 2*x11*x12*x42*xp1 - 2*x11*x12*xp1*xp2 - 2*x11*x21*x32*x42 + 2*x11*x21*x32*xp2 + 2*x11*x21*x42*xp2 - 2*x11*x21*xp2**2 - 2*x11*x22*x31*x42 + 2*x11*x22*x31*xp2 + 4*x11*x22*x32*x41 - 4*x11*x22*x32*xp1 - 4*x11*x22*x41*xp2 + 2*x11*x22*x42*xp1 + 2*x11*x22*xp1*xp2 + 2*x11*x31*x42*xp2 - 2*x11*x31*xp2**2 - 4*x11*x32*x41*xp2 + 2*x11*x32*x42*xp1 + 2*x11*x32*xp1*xp2 + 2*x11*x41*x42*xp2 + 2*x11*x41*xp2**2 - 2*x11*x42**2*xp1 - 2*x11*x42*xp1*xp2 + x12**2*x41**2 - 2*x12**2*x41*xp1 + x12**2*xp1**2 + 4*x12*x21*x31*x42 - 4*x12*x21*x31*xp2 - 2*x12*x21*x32*x41 + 2*x12*x21*x32*xp1 + 2*x12*x21*x41*xp2 - 4*x12*x21*x42*xp1 + 2*x12*x21*xp1*xp2 - 2*x12*x22*x31*x41 + 2*x12*x22*x31*xp1 + 2*x12*x22*x41*xp1 - 2*x12*x22*xp1**2 + 2*x12*x31*x41*xp2 - 4*x12*x31*x42*xp1 + 2*x12*x31*xp1*xp2 + 2*x12*x32*x41*xp1 - 2*x12*x32*xp1**2 - 2*x12*x41**2*xp2 + 2*x12*x41*x42*xp1 - 2*x12*x41*xp1*xp2 + 2*x12*x42*xp1**2 + x21**2*x32**2 - 2*x21**2*x32*xp2 + x21**2*xp2**2 - 2*x21*x22*x31*x32 + 2*x21*x22*x31*xp2 + 2*x21*x22*x32*xp1 - 2*x21*x22*xp1*xp2 + 2*x21*x31*x32*xp2 - 4*x21*x31*x42*xp2 + 2*x21*x31*xp2**2 - 2*x21*x32**2*xp1 + 2*x21*x32*x41*xp2 + 2*x21*x32*x42*xp1 - 2*x21*x32*xp1*xp2 - 2*x21*x41*xp2**2 + 2*x21*x42*xp1*xp2 + x22**2*x31**2 - 2*x22**2*x31*xp1 + x22**2*xp1**2 - 2*x22*x31**2*xp2 + 2*x22*x31*x32*xp1 + 2*x22*x31*x41*xp2 + 2*x22*x31*x42*xp1 - 2*x22*x31*xp1*xp2 - 4*x22*x32*x41*xp1 + 2*x22*x32*xp1**2 + 2*x22*x41*xp1*xp2 - 2*x22*x42*xp1**2 + x31**2*xp2**2 - 2*x31*x32*xp1*xp2 - 2*x31*x41*xp2**2 + 2*x31*x42*xp1*xp2 + x32**2*xp1**2 + 2*x32*x41*xp1*xp2 - 2*x32*x42*xp1**2 + x41**2*xp2**2 - 2*x41*x42*xp1*xp2 + x42**2*xp1**2))/(2*(x11*x22 - x11*x42 - x12*x21 + x12*x41 + x21*x32 - x22*x31 + x31*x42 - x32*x41))

  xi2b = (x11*xi1 - x11 - x21*xi1 + xp1)/(x11*xi1 - x11 - x21*xi1 - x31*xi1 + x31 + x41*xi1)
  
  return (xi1,xi2a,xi2b)

def get_xi_3d(hexahedron,p):
  [p0, p1, p2, p3, p4, p5, p6, p7] = hexahedron

  xp1 = p[0]
  xp2 = p[1]
  xp3 = p[2]
  
  x11 = p0[0]
  x12 = p0[1]
  x13 = p0[2]
  x21 = p1[0]
  x22 = p1[1]
  x23 = p1[2]
  x31 = p2[0]
  x32 = p2[1]
  x33 = p2[2]
  x41 = p3[0]
  x42 = p3[1]
  x43 = p3[2]
  x51 = p4[0]
  x52 = p4[1]
  x53 = p4[2]
  x61 = p5[0]
  x62 = p5[1]
  x63 = p5[2]
  x71 = p6[0]
  x72 = p6[1]
  x73 = p6[2]
  x81 = p7[0]
  x82 = p7[1]
  x83 = p7[2]
  
  xi3 = (2*x11*x32 - x11*x72 - x11*xp2 - 2*x12*x31 + x12*x71 + x12*xp1 + x31*x52 + x31*xp2 - x32*x51 - x32*xp1 + x51*xp2 - x52*xp1 - x71*xp2 + x72*xp1 - np.sqrt(x11**2*x72**2 - 2*x11**2*x72*xp2 + x11**2*xp2**2 - 2*x11*x12*x71*x72 + 2*x11*x12*x71*xp2 + 2*x11*x12*x72*xp1 - 2*x11*x12*xp1*xp2 - 2*x11*x31*x52*x72 + 2*x11*x31*x52*xp2 + 2*x11*x31*x72*xp2 - 2*x11*x31*xp2**2 - 2*x11*x32*x51*x72 + 2*x11*x32*x51*xp2 + 4*x11*x32*x52*x71 - 4*x11*x32*x52*xp1 - 4*x11*x32*x71*xp2 + 2*x11*x32*x72*xp1 + 2*x11*x32*xp1*xp2 + 2*x11*x51*x72*xp2 - 2*x11*x51*xp2**2 - 4*x11*x52*x71*xp2 + 2*x11*x52*x72*xp1 + 2*x11*x52*xp1*xp2 + 2*x11*x71*x72*xp2 + 2*x11*x71*xp2**2 - 2*x11*x72**2*xp1 - 2*x11*x72*xp1*xp2 + x12**2*x71**2 - 2*x12**2*x71*xp1 + x12**2*xp1**2 + 4*x12*x31*x51*x72 - 4*x12*x31*x51*xp2 - 2*x12*x31*x52*x71 + 2*x12*x31*x52*xp1 + 2*x12*x31*x71*xp2 - 4*x12*x31*x72*xp1 + 2*x12*x31*xp1*xp2 - 2*x12*x32*x51*x71 + 2*x12*x32*x51*xp1 + 2*x12*x32*x71*xp1 - 2*x12*x32*xp1**2 + 2*x12*x51*x71*xp2 - 4*x12*x51*x72*xp1 + 2*x12*x51*xp1*xp2 + 2*x12*x52*x71*xp1 - 2*x12*x52*xp1**2 - 2*x12*x71**2*xp2 + 2*x12*x71*x72*xp1 - 2*x12*x71*xp1*xp2 + 2*x12*x72*xp1**2 + x31**2*x52**2 - 2*x31**2*x52*xp2 + x31**2*xp2**2 - 2*x31*x32*x51*x52 + 2*x31*x32*x51*xp2 + 2*x31*x32*x52*xp1 - 2*x31*x32*xp1*xp2 + 2*x31*x51*x52*xp2 - 4*x31*x51*x72*xp2 + 2*x31*x51*xp2**2 - 2*x31*x52**2*xp1 + 2*x31*x52*x71*xp2 + 2*x31*x52*x72*xp1 - 2*x31*x52*xp1*xp2 - 2*x31*x71*xp2**2 + 2*x31*x72*xp1*xp2 + x32**2*x51**2 - 2*x32**2*x51*xp1 + x32**2*xp1**2 - 2*x32*x51**2*xp2 + 2*x32*x51*x52*xp1 + 2*x32*x51*x71*xp2 + 2*x32*x51*x72*xp1 - 2*x32*x51*xp1*xp2 - 4*x32*x52*x71*xp1 + 2*x32*x52*xp1**2 + 2*x32*x71*xp1*xp2 - 2*x32*x72*xp1**2 + x51**2*xp2**2 - 2*x51*x52*xp1*xp2 - 2*x51*x71*xp2**2 + 2*x51*x72*xp1*xp2 + x52**2*xp1**2 + 2*x52*x71*xp1*xp2 - 2*x52*x72*xp1**2 + x71**2*xp2**2 - 2*x71*x72*xp1*xp2 + x72**2*xp1**2))/(2*x11*x32 - 2*x11*x72 - 2*x12*x31 + 2*x12*x71 + 2*x31*x52 - 2*x32*x51 + 2*x51*x72 - 2*x52*x71)
  
  xi3b = (2*x11*x32 - x11*x72 - x11*xp2 - 2*x12*x31 + x12*x71 + x12*xp1 + x31*x52 + x31*xp2 - x32*x51 - x32*xp1 + x51*xp2 - x52*xp1 - x71*xp2 + x72*xp1 + np.sqrt(x11**2*x72**2 - 2*x11**2*x72*xp2 + x11**2*xp2**2 - 2*x11*x12*x71*x72 + 2*x11*x12*x71*xp2 + 2*x11*x12*x72*xp1 - 2*x11*x12*xp1*xp2 - 2*x11*x31*x52*x72 + 2*x11*x31*x52*xp2 + 2*x11*x31*x72*xp2 - 2*x11*x31*xp2**2 - 2*x11*x32*x51*x72 + 2*x11*x32*x51*xp2 + 4*x11*x32*x52*x71 - 4*x11*x32*x52*xp1 - 4*x11*x32*x71*xp2 + 2*x11*x32*x72*xp1 + 2*x11*x32*xp1*xp2 + 2*x11*x51*x72*xp2 - 2*x11*x51*xp2**2 - 4*x11*x52*x71*xp2 + 2*x11*x52*x72*xp1 + 2*x11*x52*xp1*xp2 + 2*x11*x71*x72*xp2 + 2*x11*x71*xp2**2 - 2*x11*x72**2*xp1 - 2*x11*x72*xp1*xp2 + x12**2*x71**2 - 2*x12**2*x71*xp1 + x12**2*xp1**2 + 4*x12*x31*x51*x72 - 4*x12*x31*x51*xp2 - 2*x12*x31*x52*x71 + 2*x12*x31*x52*xp1 + 2*x12*x31*x71*xp2 - 4*x12*x31*x72*xp1 + 2*x12*x31*xp1*xp2 - 2*x12*x32*x51*x71 + 2*x12*x32*x51*xp1 + 2*x12*x32*x71*xp1 - 2*x12*x32*xp1**2 + 2*x12*x51*x71*xp2 - 4*x12*x51*x72*xp1 + 2*x12*x51*xp1*xp2 + 2*x12*x52*x71*xp1 - 2*x12*x52*xp1**2 - 2*x12*x71**2*xp2 + 2*x12*x71*x72*xp1 - 2*x12*x71*xp1*xp2 + 2*x12*x72*xp1**2 + x31**2*x52**2 - 2*x31**2*x52*xp2 + x31**2*xp2**2 - 2*x31*x32*x51*x52 + 2*x31*x32*x51*xp2 + 2*x31*x32*x52*xp1 - 2*x31*x32*xp1*xp2 + 2*x31*x51*x52*xp2 - 4*x31*x51*x72*xp2 + 2*x31*x51*xp2**2 - 2*x31*x52**2*xp1 + 2*x31*x52*x71*xp2 + 2*x31*x52*x72*xp1 - 2*x31*x52*xp1*xp2 - 2*x31*x71*xp2**2 + 2*x31*x72*xp1*xp2 + x32**2*x51**2 - 2*x32**2*x51*xp1 + x32**2*xp1**2 - 2*x32*x51**2*xp2 + 2*x32*x51*x52*xp1 + 2*x32*x51*x71*xp2 + 2*x32*x51*x72*xp1 - 2*x32*x51*xp1*xp2 - 4*x32*x52*x71*xp1 + 2*x32*x52*xp1**2 + 2*x32*x71*xp1*xp2 - 2*x32*x72*xp1**2 + x51**2*xp2**2 - 2*x51*x52*xp1*xp2 - 2*x51*x71*xp2**2 + 2*x51*x72*xp1*xp2 + x52**2*xp1**2 + 2*x52*x71*xp1*xp2 - 2*x52*x72*xp1**2 + x71**2*xp2**2 - 2*x71*x72*xp1*xp2 + x72**2*xp1**2))/(2*x11*x32 - 2*x11*x72 - 2*x12*x31 + 2*x12*x71 + 2*x31*x52 - 2*x32*x51 + 2*x51*x72 - 2*x52*x71)
  
  xi1 = (xp1*(x13*(-xi3 + 1) + x53*xi3) - xp1*(x23*(-xi3 + 1) + x63*xi3) - xp1*(x33*(-xi3 + 1) + x73*xi3) + xp1*(x43*(-xi3 + 1) + x83*xi3) - xp2*(x11*(-xi3 + 1) + x51*xi3) + xp2*(x21*(-xi3 + 1) + x61*xi3) + xp2*(x31*(-xi3 + 1) + x71*xi3) - xp2*(x41*(-xi3 + 1) + x81*xi3) + 2*(x11*(-xi3 + 1) + x51*xi3)*(x33*(-xi3 + 1) + x73*xi3) - (x11*(-xi3 + 1) + x51*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 2*(x13*(-xi3 + 1) + x53*xi3)*(x31*(-xi3 + 1) + x71*xi3) + (x13*(-xi3 + 1) + x53*xi3)*(x41*(-xi3 + 1) + x81*xi3) - (x21*(-xi3 + 1) + x61*xi3)*(x33*(-xi3 + 1) + x73*xi3) + (x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3) + np.sqrt(xp1**2*(x13*(-xi3 + 1) + x53*xi3)**2 - 2*xp1**2*(x13*(-xi3 + 1) + x53*xi3)*(x23*(-xi3 + 1) + x63*xi3) - 2*xp1**2*(x13*(-xi3 + 1) + x53*xi3)*(x33*(-xi3 + 1) + x73*xi3) + 2*xp1**2*(x13*(-xi3 + 1) + x53*xi3)*(x43*(-xi3 + 1) + x83*xi3) + xp1**2*(x23*(-xi3 + 1) + x63*xi3)**2 + 2*xp1**2*(x23*(-xi3 + 1) + x63*xi3)*(x33*(-xi3 + 1) + x73*xi3) - 2*xp1**2*(x23*(-xi3 + 1) + x63*xi3)*(x43*(-xi3 + 1) + x83*xi3) + xp1**2*(x33*(-xi3 + 1) + x73*xi3)**2 - 2*xp1**2*(x33*(-xi3 + 1) + x73*xi3)*(x43*(-xi3 + 1) + x83*xi3) + xp1**2*(x43*(-xi3 + 1) + x83*xi3)**2 - 2*xp1*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x13*(-xi3 + 1) + x53*xi3) + 2*xp1*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x23*(-xi3 + 1) + x63*xi3) + 2*xp1*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x33*(-xi3 + 1) + x73*xi3) - 2*xp1*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 2*xp1*xp2*(x13*(-xi3 + 1) + x53*xi3)*(x21*(-xi3 + 1) + x61*xi3) + 2*xp1*xp2*(x13*(-xi3 + 1) + x53*xi3)*(x31*(-xi3 + 1) + x71*xi3) - 2*xp1*xp2*(x13*(-xi3 + 1) + x53*xi3)*(x41*(-xi3 + 1) + x81*xi3) - 2*xp1*xp2*(x21*(-xi3 + 1) + x61*xi3)*(x23*(-xi3 + 1) + x63*xi3) - 2*xp1*xp2*(x21*(-xi3 + 1) + x61*xi3)*(x33*(-xi3 + 1) + x73*xi3) + 2*xp1*xp2*(x21*(-xi3 + 1) + x61*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 2*xp1*xp2*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3) + 2*xp1*xp2*(x23*(-xi3 + 1) + x63*xi3)*(x41*(-xi3 + 1) + x81*xi3) - 2*xp1*xp2*(x31*(-xi3 + 1) + x71*xi3)*(x33*(-xi3 + 1) + x73*xi3) + 2*xp1*xp2*(x31*(-xi3 + 1) + x71*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 2*xp1*xp2*(x33*(-xi3 + 1) + x73*xi3)*(x41*(-xi3 + 1) + x81*xi3) - 2*xp1*xp2*(x41*(-xi3 + 1) + x81*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 2*xp1*(x11*(-xi3 + 1) + x51*xi3)*(x13*(-xi3 + 1) + x53*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 4*xp1*(x11*(-xi3 + 1) + x51*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x33*(-xi3 + 1) + x73*xi3) + 2*xp1*(x11*(-xi3 + 1) + x51*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 2*xp1*(x11*(-xi3 + 1) + x51*xi3)*(x33*(-xi3 + 1) + x73*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 2*xp1*(x11*(-xi3 + 1) + x51*xi3)*(x43*(-xi3 + 1) + x83*xi3)**2 - 2*xp1*(x13*(-xi3 + 1) + x53*xi3)**2*(x41*(-xi3 + 1) + x81*xi3) + 2*xp1*(x13*(-xi3 + 1) + x53*xi3)*(x21*(-xi3 + 1) + x61*xi3)*(x33*(-xi3 + 1) + x73*xi3) - 4*xp1*(x13*(-xi3 + 1) + x53*xi3)*(x21*(-xi3 + 1) + x61*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 2*xp1*(x13*(-xi3 + 1) + x53*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3) + 2*xp1*(x13*(-xi3 + 1) + x53*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x41*(-xi3 + 1) + x81*xi3) - 4*xp1*(x13*(-xi3 + 1) + x53*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 2*xp1*(x13*(-xi3 + 1) + x53*xi3)*(x33*(-xi3 + 1) + x73*xi3)*(x41*(-xi3 + 1) + x81*xi3) + 2*xp1*(x13*(-xi3 + 1) + x53*xi3)*(x41*(-xi3 + 1) + x81*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 2*xp1*(x21*(-xi3 + 1) + x61*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x33*(-xi3 + 1) + x73*xi3) - 2*xp1*(x21*(-xi3 + 1) + x61*xi3)*(x33*(-xi3 + 1) + x73*xi3)**2 + 2*xp1*(x21*(-xi3 + 1) + x61*xi3)*(x33*(-xi3 + 1) + x73*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 2*xp1*(x23*(-xi3 + 1) + x63*xi3)**2*(x31*(-xi3 + 1) + x71*xi3) + 2*xp1*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x33*(-xi3 + 1) + x73*xi3) + 2*xp1*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 4*xp1*(x23*(-xi3 + 1) + x63*xi3)*(x33*(-xi3 + 1) + x73*xi3)*(x41*(-xi3 + 1) + x81*xi3) + xp2**2*(x11*(-xi3 + 1) + x51*xi3)**2 - 2*xp2**2*(x11*(-xi3 + 1) + x51*xi3)*(x21*(-xi3 + 1) + x61*xi3) - 2*xp2**2*(x11*(-xi3 + 1) + x51*xi3)*(x31*(-xi3 + 1) + x71*xi3) + 2*xp2**2*(x11*(-xi3 + 1) + x51*xi3)*(x41*(-xi3 + 1) + x81*xi3) + xp2**2*(x21*(-xi3 + 1) + x61*xi3)**2 + 2*xp2**2*(x21*(-xi3 + 1) + x61*xi3)*(x31*(-xi3 + 1) + x71*xi3) - 2*xp2**2*(x21*(-xi3 + 1) + x61*xi3)*(x41*(-xi3 + 1) + x81*xi3) + xp2**2*(x31*(-xi3 + 1) + x71*xi3)**2 - 2*xp2**2*(x31*(-xi3 + 1) + x71*xi3)*(x41*(-xi3 + 1) + x81*xi3) + xp2**2*(x41*(-xi3 + 1) + x81*xi3)**2 - 2*xp2*(x11*(-xi3 + 1) + x51*xi3)**2*(x43*(-xi3 + 1) + x83*xi3) + 2*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x13*(-xi3 + 1) + x53*xi3)*(x41*(-xi3 + 1) + x81*xi3) + 2*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x21*(-xi3 + 1) + x61*xi3)*(x33*(-xi3 + 1) + x73*xi3) + 2*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x21*(-xi3 + 1) + x61*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 2*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3) - 4*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x41*(-xi3 + 1) + x81*xi3) + 2*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 4*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x33*(-xi3 + 1) + x73*xi3)*(x41*(-xi3 + 1) + x81*xi3) + 2*xp2*(x11*(-xi3 + 1) + x51*xi3)*(x41*(-xi3 + 1) + x81*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 4*xp2*(x13*(-xi3 + 1) + x53*xi3)*(x21*(-xi3 + 1) + x61*xi3)*(x31*(-xi3 + 1) + x71*xi3) + 2*xp2*(x13*(-xi3 + 1) + x53*xi3)*(x21*(-xi3 + 1) + x61*xi3)*(x41*(-xi3 + 1) + x81*xi3) + 2*xp2*(x13*(-xi3 + 1) + x53*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x41*(-xi3 + 1) + x81*xi3) - 2*xp2*(x13*(-xi3 + 1) + x53*xi3)*(x41*(-xi3 + 1) + x81*xi3)**2 - 2*xp2*(x21*(-xi3 + 1) + x61*xi3)**2*(x33*(-xi3 + 1) + x73*xi3) + 2*xp2*(x21*(-xi3 + 1) + x61*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3) + 2*xp2*(x21*(-xi3 + 1) + x61*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x33*(-xi3 + 1) + x73*xi3) - 4*xp2*(x21*(-xi3 + 1) + x61*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 2*xp2*(x21*(-xi3 + 1) + x61*xi3)*(x33*(-xi3 + 1) + x73*xi3)*(x41*(-xi3 + 1) + x81*xi3) - 2*xp2*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3)**2 + 2*xp2*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x41*(-xi3 + 1) + x81*xi3) + (x11*(-xi3 + 1) + x51*xi3)**2*(x43*(-xi3 + 1) + x83*xi3)**2 - 2*(x11*(-xi3 + 1) + x51*xi3)*(x13*(-xi3 + 1) + x53*xi3)*(x41*(-xi3 + 1) + x81*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 2*(x11*(-xi3 + 1) + x51*xi3)*(x21*(-xi3 + 1) + x61*xi3)*(x33*(-xi3 + 1) + x73*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 2*(x11*(-xi3 + 1) + x51*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 4*(x11*(-xi3 + 1) + x51*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x33*(-xi3 + 1) + x73*xi3)*(x41*(-xi3 + 1) + x81*xi3) + (x13*(-xi3 + 1) + x53*xi3)**2*(x41*(-xi3 + 1) + x81*xi3)**2 + 4*(x13*(-xi3 + 1) + x53*xi3)*(x21*(-xi3 + 1) + x61*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 2*(x13*(-xi3 + 1) + x53*xi3)*(x21*(-xi3 + 1) + x61*xi3)*(x33*(-xi3 + 1) + x73*xi3)*(x41*(-xi3 + 1) + x81*xi3) - 2*(x13*(-xi3 + 1) + x53*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x41*(-xi3 + 1) + x81*xi3) + (x21*(-xi3 + 1) + x61*xi3)**2*(x33*(-xi3 + 1) + x73*xi3)**2 - 2*(x21*(-xi3 + 1) + x61*xi3)*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3)*(x33*(-xi3 + 1) + x73*xi3) + (x23*(-xi3 + 1) + x63*xi3)**2*(x31*(-xi3 + 1) + x71*xi3)**2))/(2*(x11*(-xi3 + 1) + x51*xi3)*(x33*(-xi3 + 1) + x73*xi3) - 2*(x11*(-xi3 + 1) + x51*xi3)*(x43*(-xi3 + 1) + x83*xi3) - 2*(x13*(-xi3 + 1) + x53*xi3)*(x31*(-xi3 + 1) + x71*xi3) + 2*(x13*(-xi3 + 1) + x53*xi3)*(x41*(-xi3 + 1) + x81*xi3) - 2*(x21*(-xi3 + 1) + x61*xi3)*(x33*(-xi3 + 1) + x73*xi3) + 2*(x21*(-xi3 + 1) + x61*xi3)*(x43*(-xi3 + 1) + x83*xi3) + 2*(x23*(-xi3 + 1) + x63*xi3)*(x31*(-xi3 + 1) + x71*xi3) - 2*(x23*(-xi3 + 1) + x63*xi3)*(x41*(-xi3 + 1) + x81*xi3))

  xi2 = (x11*xi1*xi3 - x11*xi1 - x11*xi3 + x11 - x21*xi1*xi3 + x21*xi1 - x51*xi1*xi3 + x51*xi3 + x61*xi1*xi3 - xp1)/(x11*xi1*xi3 - x11*xi1 - x11*xi3 + x11 - x21*xi1*xi3 + x21*xi1 - x31*xi1*xi3 + x31*xi1 + x31*xi3 - x31 + x41*xi1*xi3 - x41*xi1 - x51*xi1*xi3 + x51*xi3 + x61*xi1*xi3 + x71*xi1*xi3 - x71*xi3 - x81*xi1*xi3)
  return (xi1,xi2,xi3,xi3b)

np.random.seed(1)
max_factor = 0
def point_is_in_tetrahedron(tetrahedron,correct_orientation,p):
  [p3, p0, p1, p2] = tetrahedron
  global max_factor
  
  debug = False
  
  xp1 = p[0]
  xp2 = p[1]
  xp3 = p[2]
  
  x11 = p0[0]
  x12 = p0[1]
  x13 = p0[2]
  x21 = p1[0]
  x22 = p1[1]
  x23 = p1[2]
  x31 = p2[0]
  x32 = p2[1]
  x33 = p2[2]
  x41 = p3[0]
  x42 = p3[1]
  x43 = p3[2]
  
  if debug:
    print ""
    print "----"
    print "tetrahedron:", tetrahedron, ", p:",p
  det = (x11 - x41)*(x22 - x42)*(x33 - x43) - (x11 - x41)*(x23 - x43)*(x32 - x42) - (x12 - x42)*(x21 - x41)*(x33 - x43) + (x12 - x42)*(x23 - x43)*(x31 - x41) + (x13 - x43)*(x21 - x41)*(x32 - x42) - (x13 - x43)*(x22 - x42)*(x31 - x41)
  xi1 = 1/det * ((-x41 + xp1)*((x22 - x42)*(x33 - x43) - (x23 - x43)*(x32 - x42)) + (-x42 + xp2)*(-(x21 - x41)*(x33 - x43) + (x23 - x43)*(x31 - x41)) + (-x43 + xp3)*((x21 - x41)*(x32 - x42) - (x22 - x42)*(x31 - x41)))
  xi2 = 1/det * ((-x41 + xp1)*(-(x12 - x42)*(x33 - x43) + (x13 - x43)*(x32 - x42)) + (-x42 + xp2)*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) + (-x43 + xp3)*(-(x11 - x41)*(x32 - x42) + (x12 - x42)*(x31 - x41)))
  xi3 = 1/det * ((-x41 + xp1)*((x12 - x42)*(x23 - x43) - (x13 - x43)*(x22 - x42)) + (-x42 + xp2)*(-(x11 - x41)*(x23 - x43) + (x13 - x43)*(x21 - x41)) + (-x43 + xp3)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41)))

  if debug:
    print "det: ", det, ",xi1: ",xi1,", xi2:",xi2, ",xi3:",xi3
    
  factor = 3-xi1**2-xi2**2-xi3**2
  #factor = (np.sqrt(3)-np.linalg.norm(np.array([xi1,xi2,xi3])))
 
  point_is_inside = (xi1 >= 0 and xi2 >= 0 and xi3 >= 0)
  max_factor = max(factor,max_factor)
  if not correct_orientation[0]:
    xi1 = 1. - xi1
  if not correct_orientation[1]:
    xi2 = 1. - xi2
  if not correct_orientation[2]:
    xi3 = 1. - xi3
    
  if debug:
    print "t matrix: "
    tmat = np.array([[x11-x41,x21-x41,x31-x41],[x12-x42,x22-x42,x32-x42],[x13-x43,x23-x43,x33-x43]])
    print tmat
    print "adj: "
    adj = np.array([[(x22 - x42)*(x33 - x43) - (x23 - x43)*(x32 - x42), -(x21 - x41)*(x33 - x43) + (x23 - x43)*(x31 - x41), (x21 - x41)*(x32 - x42) - (x22 - x42)*(x31 - x41)], [-(x12 - x42)*(x33 - x43) + (x13 - x43)*(x32 - x42), (x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41), -(x11 - x41)*(x32 - x42) + (x12 - x42)*(x31 - x41)], [(x12 - x42)*(x23 - x43) - (x13 - x43)*(x22 - x42), -(x11 - x41)*(x23 - x43) + (x13 - x43)*(x21 - x41), (x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41)]])
    print adj
    print "det: ", det, ",xi1: ",xi1,", xi2:",xi2, ",xi3:",xi3,",   factor:",factor
    print ""
    print "check:"
    xi = np.array([[xi1],[xi2],[xi3]])
    xi = np.array([[1.0],[0.0],[0.0]])
    print tmat.dot(xi)
    print ""
    print tmat.dot(xi),"=",p-p3
    print "p-p3:",(p - p3)
    tinv = np.array([[(-(-(-x12 + x42)*((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41)) + (-x13 + x43)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41)))*(-(x21 - x41)*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)) + (x31 - x41)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))) + (((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)))*((x11 - x41)*(x22 - x42) - (-x12 + x42)*(x21 - x41) - (x12 - x42)*(x21 - x41)))/((x11 - x41)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)))), (-(x11 - x41)*(x21 - x41)*(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41))) + (x11 - x41)*((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*(-(x21 - x41)*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)) + (x31 - x41)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))))/((x11 - x41)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)))), -(-(x21 - x41)*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)) + (x31 - x41)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41)))/(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)))], [((-x12 + x42)*(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41))) - ((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41))*(-(-x12 + x42)*((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41)) + (-x13 + x43)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))))/(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)))), ((x11 - x41)*((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)) + (x11 - x41)*(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41))))/(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)))), -(x11 - x41)*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41))/(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)))], [(-(-x12 + x42)*((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41)) + (-x13 + x43)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41)))/(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41))), -(x11 - x41)*((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))/(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41))), (x11 - x41)*((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))/(((x11 - x41)*(x22 - x42) - (x12 - x42)*(x21 - x41))*((x11 - x41)*(x33 - x43) - (x13 - x43)*(x31 - x41)) - ((x11 - x41)*(x23 - x43) - (x13 - x43)*(x21 - x41))*((x11 - x41)*(x32 - x42) - (x12 - x42)*(x31 - x41)))]])
    print "t^-1:"
    print tinv
    det = (x11 - x41)*(x22 - x42)*(x33 - x43) - (x11 - x41)*(x23 - x43)*(x32 - x42) - (x12 - x42)*(x21 - x41)*(x33 - x43) + (x12 - x42)*(x23 - x43)*(x31 - x41) + (x13 - x43)*(x21 - x41)*(x32 - x42) - (x13 - x43)*(x22 - x42)*(x31 - x41)
    print "det:",det
    print "adj/det:"
    print adj/det
    
    xir = tinv.dot(np.reshape(p-p3,(3,1)))
    print "t^-1 (p-p3)="
    print xir
    print "="
    print (adj/det).dot(np.reshape(p-p3,(3,1)))
    
  return (point_is_inside, (xi1, xi2, xi3), factor)
  
def get_xi_3d2(hexahedron,p):
  [p0, p1, p2, p3, p4, p5, p6, p7] = hexahedron

  debug = False

  xi_sum = np.zeros(3)
  denom = 0
  no_factor = True
  only_some = False

  # p0
  r = point_is_in_tetrahedron([p0, p1, p2, p4], [True,True,True], p)
  xi = np.array(r[1])
  if not r[0]:
    if debug:
      print "p0 out"
  if debug:
    print "0 xi: ",xi
    
  factor = 1./(0.1+np.linalg.norm(p-p0))
  if no_factor:
    factor = 1.0
  xi_sum += xi*factor
  denom += factor
    
  # p1
  if not only_some:
    r = point_is_in_tetrahedron([p1, p0, p5, p3], [False,True,True], p)
    xi = np.array(r[1])
    if not r[0]:
      if debug:
        print "p1 out"
    xi[1],xi[2] = xi[2],xi[1]
    if debug:
      print "1 xi: ",xi
      
    factor = 1./(0.1+np.linalg.norm(p-p1))
    if no_factor:
      factor = 1.0
    xi_sum += xi*factor
    denom += factor
      
    # p2
    r = point_is_in_tetrahedron([p2, p3, p6, p0], [True,True,False], p)
    xi = np.array(r[1])
    if not r[0]:
      if debug:
        print "p2 out"
    xi[1],xi[2] = xi[2],xi[1]
    if debug:
      print "2 xi: ",xi
      
    factor = 1./(0.1+np.linalg.norm(p-p2))
    if no_factor:
      factor = 1.0
    xi_sum += xi*factor
    denom += factor
    
  # p3
  r = point_is_in_tetrahedron([p3, p2, p1, p7], [False,False,True], p)
  xi = np.array(r[1])
  if not r[0]:
    if debug:
      print "p3 out"
  if debug:
    print "3 xi: ",xi
  
  factor = 1./(0.1+np.linalg.norm(p-p3))
  if no_factor:
    factor = 1.0
  xi_sum += xi*factor
  denom += factor
  
  # p4
  if not only_some:
    r = point_is_in_tetrahedron([p4, p5, p0, p6], [True,False,True], p)
    xi = np.array(r[1])
    if not r[0]:
      if debug:
        print "p4 out"
    xi[1],xi[2] = xi[2],xi[1]
    if debug:
      print "4 xi: ",xi
    
    factor = 1./(0.1+np.linalg.norm(p-p4))
    if no_factor:
      factor = 1.0
    xi_sum += xi*factor
    denom += factor
  
  # p5
  r = point_is_in_tetrahedron([p5, p4, p7, p1], [False,True,False], p)
  xi = np.array(r[1])
  if not r[0]:
    if debug:
      print "p5 out"
  if debug:
    print "5 xi: ",xi
  
  factor = 1./(0.1+np.linalg.norm(p-p5))
  if no_factor:
    factor = 1.0
  xi_sum += xi*factor
  denom += factor
  
  # p6
  r = point_is_in_tetrahedron([p6, p7, p4, p2], [True,False,False], p)
  xi = np.array(r[1])
  if not r[0]:
    if debug:
      print "p6 out"
  if debug:
    print "6 xi: ",xi

  factor = 1./(0.1+np.linalg.norm(p-p6))
  if no_factor:
    factor = 1.0
  xi_sum += xi*factor
  denom += factor
  
  # p7
  if not only_some:
    r = point_is_in_tetrahedron([p7, p6, p3, p5], [False,False,False], p)
    xi = np.array(r[1])
    if not r[0]:
      if debug:
        print "p7 out"
    xi[1],xi[2] = xi[2],xi[1]
    if debug:
      print "7 xi: ",xi

    factor = 1./(0.1+np.linalg.norm(p-p7))
    if no_factor:
      factor = 1.0
    xi_sum += xi*factor
    denom += factor
    
  xi_sum /= denom 
  
  if debug:
    print "final xi: ", xi_sum
  
  eps = 1e-12
  if (0.0-eps <= xi_sum[0] <= 1.0+eps) and (0.0-eps <= xi_sum[1] <= 1.0+eps) and (0.0-eps <= xi_sum[2] <= 1.0+eps):
    if debug:
      print "inside"
    return xi_sum
  else:
    if debug:
      print "outside"
    return None

def point_is_in_element(hexahedron,p):
  [p0, p1, p2, p3, p4, p5, p6, p7] = hexahedron
  
  # bottom [p0,p1,p3,p2]
  a30 = (-p3+p0)
  a01 = (-p0+p1)
  a12 = (-p1+p2)
  a32 = (-p3+p2)
  a20 = (-p2+p0)
  
  v0 = np.cross(a30, a01).dot(-p0+p) >= 0
  v1 = np.cross(a01, a12).dot(-p1+p) >= 0
  v2 = np.cross(a12, a32).dot(-p3+p) >= 0
  v3 = np.cross(a32, a20).dot(-p2+p) >= 0
  
  # top [p4,p6,p7,p5]
  a74 = (-p7+p4)
  a46 = (-p4+p6)
  a65 = (-p6+p5)
  a75 = (-p7+p5)
  a54 = (-p5+p4)
  
  v4 = np.cross(a74, a46).dot(-p4+p) >= 0
  v5 = np.cross(a46, a65).dot(-p6+p) >= 0
  v6 = np.cross(a65, a75).dot(-p7+p) >= 0
  v7 = np.cross(a75, a54).dot(-p5+p) >= 0
  
  # right [p1,p5,p7,p3]
  a71 = (-p7+p1)
  a15 = (-p1+p5)
  a53 = (-p5+p3)
  a73 = (-p7+p3)
  a31 = (-p3+p1)
  
  v8 = np.cross(a71, a15).dot(-p1+p) >= 0
  v9 = np.cross(a15, a53).dot(-p5+p) >= 0
  v10 = np.cross(a53, a73).dot(-p7+p) >= 0
  v11 = np.cross(a73, a31).dot(-p3+p) >= 0
  
  # left [p0,p2,p6,p4]
  a60 = (-p6+p0)
  a02 = (-p0+p2)
  a24 = (-p2+p4)
  a64 = (-p6+p4)
  a40 = (-p4+p0)
  
  v12 = np.cross(a60, a02).dot(-p0+p) >= 0
  v13 = np.cross(a02, a24).dot(-p2+p) >= 0
  v14 = np.cross(a24, a64).dot(-p6+p) >= 0
  v15 = np.cross(a64, a40).dot(-p4+p) >= 0
  
  # front [p0,p4,p5,p1]
  a50 = (-p5+p0)
  a04 = (-p0+p4)
  a41 = (-p4+p1)
  a51 = (-p5+p1)
  a10 = (-p1+p0)
  
  v16 = np.cross(a50, a04).dot(-p0+p) >= 0
  v17 = np.cross(a04, a41).dot(-p4+p) >= 0
  v18 = np.cross(a41, a51).dot(-p5+p) >= 0
  v19 = np.cross(a51, a10).dot(-p1+p) >= 0
  
  # back [p2,p3,p7,p6]
  a72 = (-p7+p2)
  a23 = (-p2+p3)
  a36 = (-p3+p6)
  a76 = (-p7+p6)
  a62 = (-p6+p2)
  
  v20 = np.cross(a72, a23).dot(-p2+p) >= 0
  v21 = np.cross(a23, a36).dot(-p3+p) >= 0
  v22 = np.cross(a36, a76).dot(-p7+p) >= 0
  v23 = np.cross(a76, a62).dot(-p6+p) >= 0
  
  b1 = point_is_in_front_of_quadrilateral([p0,p1,p3,p2],p)   # bottom
  b2 = point_is_in_front_of_quadrilateral([p4,p6,p7,p5],p)  # top
  b3 = point_is_in_front_of_quadrilateral([p1,p5,p7,p3],p)   # right
  b4 = point_is_in_front_of_quadrilateral([p0,p2,p6,p4],p)   # left
  b5 = point_is_in_front_of_quadrilateral([p0,p4,p5,p1],p)   # front
  b6 = point_is_in_front_of_quadrilateral([p2,p3,p7,p6],p)   # back
  
  is_inside0 = v0 and v1 and v2 and v3 and v4 and v5 and v6 and v7 and v8 and v9 and v10 and v11 and v12 and v13 and v14 and v15 and v16 and v17 and v18 and v19 and v20 and v21 and v22 and v23
  is_inside = b1 and b2 and b3 and b4 and b5 and b6
  
  if is_inside0 != is_inside:
    print "b: ",b1,b2,b3,b4,b5,b6
    print "v: ",v0, v1, v2, v3, ",", v4, v5, v6, v7, ",", v8, v9, v10, v11, ",", v12, v13, v14, v15, ",",  v16, v17, v18, v19, ",", v20, v21, v22, v23
    print "error!"
  
  debug = True
  if debug:
    print "v: ",v0, v1, v2, v3, ",", v4, v5, v6, v7, ",", v8, v9, v10, v11, ",", v12, v13, v14, v15, ",",  v16, v17, v18, v19  , ",", v20, v21, v22, v23
    
    print ""
    print "point ",p
    if b1 and b2 and b3 and b4 and b5 and b6:
      print "inside"
    else:
      print "outside"
      
    pp = p
    p = [p0, p1, p2, p3, p4, p5, p6, p7]

    import stl
    from stl import mesh

    out_3d_mesh_triangles = [
    [p[1],p[0],p[2]],[p[1],p[2],p[3]],  # bottom
    [p[0],p[3],p[1]],[p[0],p[2],p[3]],  # bottom
    [p[4],p[5],p[7]],[p[4],p[7],p[6]],  # top
    [p[0],p[1],p[5]],[p[0],p[5],p[4]],  # front
    [p[2],p[7],p[3]],[p[2],p[6],p[7]],  # back
    [p[2],p[0],p[4]],[p[2],p[4],p[6]],  # left
    [p[1],p[3],p[7]],[p[1],p[7],p[5]],  # right
    [pp, pp+np.array([0.1,0.0,0.0]), pp-np.array([0.1,0.0,0.0])],
    [pp, pp+np.array([0.0,0.1,0.0]), pp-np.array([0.0,0.1,0.0])],
    [pp, pp+np.array([0.0,0.0,0.1]), pp-np.array([0.0,0.0,0.1])],
    ]

    # write debugging output stl meshes
    def write_stl(triangles, outfile, description):
      # create output mesh
      n_triangles = len(triangles)

      # Create the mesh
      out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
      for i, f in enumerate(triangles):
        out_mesh.vectors[i] = f
        #for j in range(3):
          #print "set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
          
      #out_mesh.update_normals()

      out_mesh.save(outfile, mode=stl.Mode.ASCII)
      print "saved {} triangles to \"{}\" ({})".format(n_triangles,outfile,description)

    write_stl(out_3d_mesh_triangles,   "auat.stl", "aut")

  
  
  return is_inside
  
def point_is_in_front_of_quadrilateral(quadrilateral,p):
  [p0, p1, p2, p3] = quadrilateral
  
  v0 = np.cross((-p2+p0), (-p0+p1)).dot(-p0+p)
  v1 = np.cross((-p0+p1), (-p1+p3)).dot(-p1+p)
  v2 = np.cross((-p1+p2), (-p2+p3)).dot(-p2+p)
  v3 = np.cross((-p2+p3), (-p3+p0)).dot(-p3+p)
  
  c2 = np.cross((-p1+p2), (-p2+p3))
  c2 = c2/np.linalg.norm(c2)
  c22 = -p2+p
  c22 = c22/np.linalg.norm(c22)
  
  
  v2 = c2.dot(c22)
  print "angle=",np.arccos(v2)*180./np.pi
  
  
  debug = True
  if debug:
    print "quad ",quadrilateral
    print "c2=",c2,",p22=",c22
    print "v2 = ",(-p1+p2),"x",(-p2+p3)," (=",np.cross((-p1+p2), (-p2+p3)),"), dot",-p2+p
    print ", v0:",v0,", v1:",v1,", v2:",v2,", v3:",v3
  
  return v0 >= 0 and v1 >= 0 and v2 >= 0 and v3 >= 0

if True:
#point (58.9434,146.219,37), element 1
# p0 ((61.6746,146.275,37),
# p1 (60.3037,146.374,37),
# p2 (62.1413,145.011,37),
# p3 (60.6606,145.111,37),
# p4 (71.9646,148.898,50.8421),
# p5 (69.9016,149.867,50.8421),
# p6 (71.7574,146.739,50.8421),
# p7 (69.8146,147.667,50.8421))
#DEBUG:   0 xi: (2.06187,0.204832,4.08175e-16)
#DEBUG:   1 xi: (2.04524,4.08645e-16,0.203918)
#DEBUG:   2 xi: (1.90644,4.12565e-16,0.196281)
#DEBUG:   3 xi: (1.96677,0.199601,4.10861e-16)
#DEBUG:   4 xi: (1.26432,0,0.592938)
#DEBUG:   5 xi: (1.64441,0.353881,3.33067e-16)
#DEBUG:   6 xi: (1.63102,1.14093,5.55112e-16)
#DEBUG:   7 xi: (1.88966,5.55112e-16,0.871531)

#VERB3: pointIsInElement, point (63.0191,146.036,37), element 0((63.0191,146.036,37),(61.6746,146.275,37),(64.2323,145.398,37),(62.1413,145.011,37),(73.7444,147.455,50.8421),(71.9646,148.898,50.8421),(75.2802,145.762,50.8421),(71.7574,146.739,50.8421))
#VERB3:    xi: (0,0,0)
#VERB3:    xi: (1.11022e-16,-7.1741e-32,-3.73501e-18)
#VERB3:    xi: (0,0,0)
#VERB3:    xi: (0.43771,0.361389,3.27812e-16)
#VERB3:    xi: (0,0,0)
#VERB3:    xi: (0.287173,-0.365592,7.77156e-16)
#VERB3: p5 out


  # test
  #DEBUG: pointIsInElement, point 
  #  p (59.368,144.955,37), element 7
  # p0 (59.368,144.955,37),
  # p1 (57.0694,144.573,37),
  # p2 (59.3694,144.08,37),
  # p3 (56.4743,143.333,37),
  # p4 (67.8127,148.044,50.8421),
  # p5 (64.6207,148.895,50.8421),
  # p6 (67.3157,146.434,50.8421),
  # p7 (63.9783,146.701,50.8421)
  #DEBUG:    xi: (0,0,0)
  #DEBUG:    xi: (3.06205,0,-11.8288)
  #DEBUG: p1 out
  # [p1, p0, p5, p3]
  #VERB3:  isInside: 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1

# ((),(),(),(),(),(),(),()

  p0 = np.array((72.253,151.97,120.053))
  p1 = np.array((57.9479,156.418,120.053))
  p2 = np.array((71.5834,139.519,120.053))
  p3 = np.array((64.2421,147.55,120.053))
  p4 = np.array((72.9941,155.63,133.895))
  p5 = np.array((57.1602,159.617,133.895))
  p6 = np.array((71.1095,141.566,133.895))
  p7 = np.array((63.6357,150.208,133.895))
  pp = np.array([70,150,122])
  
  p = [p0, p1, p2, p3, p4, p5, p6, p7]

  import stl
  from stl import mesh

  out_3d_mesh_triangles = [
  [p[0],p[3],p[1]],[p[0],p[2],p[3]],  # bottom
  [p[4],p[5],p[7]],[p[4],p[7],p[6]],  # top
  [p[0],p[1],p[5]],[p[0],p[5],p[4]],  # front
  [p[2],p[7],p[3]],[p[2],p[6],p[7]],  # back
  [p[2],p[0],p[4]],[p[2],p[4],p[6]],  # left
  [p[1],p[3],p[7]],[p[1],p[7],p[5]],  # right
  [pp, pp+np.array([0.1,0.0,0.0]), pp+np.array([0.2,0.0,0.0])],
  ]

  # write debugging output stl meshes
  def write_stl(triangles, outfile, description):
    # create output mesh
    n_triangles = len(triangles)

    # Create the mesh
    out_mesh = mesh.Mesh(np.zeros(n_triangles, dtype=mesh.Mesh.dtype))
    for i, f in enumerate(triangles):
      out_mesh.vectors[i] = f
      #for j in range(3):
        #print "set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
        
    #out_mesh.update_normals()

    out_mesh.save(outfile, mode=stl.Mode.ASCII)
    print "saved {} triangles to \"{}\" ({})".format(n_triangles,outfile,description)

  write_stl(out_3d_mesh_triangles,   "aut.stl", "aut")

  p  = pp


  point_is_in_element([p0, p1, p2, p3, p4, p5, p6, p7], p)
    
  sys.exit(0)

# test
factor = 0.050
#factor = 0
p0 = np.array([0.0, 0.0, 0.0]) + np.random.rand(3)*factor
p1 = np.array([1.0, 0.0, 0.0]) + np.random.rand(3)*factor
p2 = np.array([0.0, 1.0, 0.0]) + np.random.rand(3)*factor
p3 = np.array([1.0, 1.0, 0.0]) + np.random.rand(3)*factor
p4 = np.array([0.0, 0.0, 1.0]) + np.random.rand(3)*factor
p5 = np.array([1.0, 0.0, 1.0]) + np.random.rand(3)*factor
p6 = np.array([0.0, 1.0, 1.0]) + np.random.rand(3)*factor
p7 = np.array([1.0, 1.0, 1.0]) + np.random.rand(3)*factor
p  = np.array([0.3, 0.4, 0.2])

print "3d"


xis = np.random.rand(100000,3)
#xis = 

error_sum = 0
n = 0
for xi in xis:
  (xi1,xi2,xi3) = xi
  p = (1-xi1)*(1-xi2)*(1-xi3)*p0 + xi1*(1-xi2)*(1-xi3)*p1 + (1-xi1)*xi2*(1-xi3)*p2 + xi1*xi2*(1-xi3)*p3 + (1-xi1)*(1-xi2)*xi3*p4 + xi1*(1-xi2)*xi3*p5 + (1-xi1)*xi2*xi3*p6 + xi1*xi2*xi3*p7
  
  xi_comp = get_xi_3d2([p0, p1, p2, p3, p4, p5, p6, p7], p)
  
  in_el0 = xi_comp is not None
  in_el1 = point_is_in_element([p0, p1, p2, p3, p4, p5, p6, p7], p)
  
  if in_el0 != in_el1:
    print "break"
    #break
  
  if xi_comp is None:
    #print "None"
    continue
  
  error = np.linalg.norm(xi_comp-xi)
  error_sum += error
  n += 1
  #if error < 1e-12:
  #  print xi,"ok"
  #else:
  #  print xi,"failed, error: ", error,", computed: ",xi_comp,", point: ",p

print "avg error: ", error_sum/n
print "max_factor:",max_factor
#print "2d"
#print get_xi_2d([p0, p1, p2, p3], p)

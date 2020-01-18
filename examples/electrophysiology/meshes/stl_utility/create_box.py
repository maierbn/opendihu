#!/usr/bin/env python

import numpy as np
from numpy import sin, cos
import random
import sys
import csv
import struct

# read in arguments
if len(sys.argv) < 5:
  print "usage: {} <filename> <x> <y> <z>\n Transform vertices such that center of mass is at [x,y,z].".format(sys.argv[0])
  sys.exit(0)

filename_out = sys.argv[1]
x = 0.0
y = 0.0
z = 0.0

if len(sys.argv) == 5:
 x = float(sys.argv[2])
 y = float(sys.argv[3])
 z = float(sys.argv[4])
       
print("\nOutput file: {}, dimensions: [{},{},{}]".format(filename_out, x, y, z))

def write_triangle(fout, p0, p1, p2):
  
  #print "p0: "+str(p0[0])+" "+str(p0[1])+" "+str(p0[2])
  #print "p1: "+str(p1[0])+" "+str(p1[1])+" "+str(p1[2])
  #print "p2: "+str(p2[0])+" "+str(p2[1])+" "+str(p2[2])
  
  # compute normal
  normal = np.cross(p1-p0, p2-p0)
  if abs(np.linalg.norm(normal)) < 1e-14:
    return
  normal = normal / np.linalg.norm(normal)
  
  scale = 1    # scale from meters to millimeters

  fout.write("facet normal "+str(normal[0])+" "+str(normal[1])+" "+str(normal[2])+"\n")
  fout.write("\touter loop\n")
  fout.write("\t\tvertex "+str(p0[0]*scale)+" "+str(p0[1]*scale)+" "+str(p0[2]*scale)+"\n")
  fout.write("\t\tvertex "+str(p1[0]*scale)+" "+str(p1[1]*scale)+" "+str(p1[2]*scale)+"\n")
  fout.write("\t\tvertex "+str(p2[0]*scale)+" "+str(p2[1]*scale)+" "+str(p2[2]*scale)+"\n")
  fout.write("\tendloop\n")
  fout.write("endfacet\n")
  
"""   p2 p3 
      p0 p1 """
def write_rectangle(fout, p0, p1, p2, p3):
  write_triangle(fout, p0, p1, p3)
  write_triangle(fout, p0, p3, p2)
  
def write_cuboid(fout, xa, xb, ya, yb, za, zb):
  
  # front
  p0 = np.array([xa, ya, za])
  p1 = np.array([xb, ya, za])
  p2 = np.array([xa, ya, zb])
  p3 = np.array([xb, ya, zb])
  write_rectangle(fout, p0, p1, p2, p3)
    
  # back
  p0 = np.array([xb, yb, za])
  p1 = np.array([xa, yb, za])
  p2 = np.array([xb, yb, zb])
  p3 = np.array([xa, yb, zb])
  write_rectangle(fout, p0, p1, p2, p3)
  
  # left
  p0 = np.array([xa, yb, za])
  p1 = np.array([xa, ya, za])
  p2 = np.array([xa, yb, zb])
  p3 = np.array([xa, ya, zb])
  write_rectangle(fout, p0, p1, p2, p3)
  
  # right
  p0 = np.array([xb, ya, za])
  p1 = np.array([xb, yb, za])
  p2 = np.array([xb, ya, zb])
  p3 = np.array([xb, yb, zb])
  write_rectangle(fout, p0, p1, p2, p3)
    
  # bottom
  p0 = np.array([xa, yb, za])
  p1 = np.array([xb, yb, za])
  p2 = np.array([xa, ya, za])
  p3 = np.array([xb, ya, za])
  write_rectangle(fout, p0, p1, p2, p3)
  
  # top
  p0 = np.array([xb, yb, zb])
  p1 = np.array([xa, yb, zb])
  p2 = np.array([xb, ya, zb])
  p3 = np.array([xa, ya, zb])
  write_rectangle(fout, p0, p1, p2, p3)

# create file
fout = open(filename_out, 'w')
fout.write("solid model\n")

xa = 0
ya = 0
za = 0
xb = x
yb = y
zb = z
write_cuboid(fout, xa, xb, ya, yb, za, zb)

fout.close()

print("File {} written.".format(filename_out))

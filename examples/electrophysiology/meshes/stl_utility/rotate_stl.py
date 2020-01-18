#!/usr/bin/env python

# @file rotate_stl.py
# 
# (c) Robert Bosch GmbH 2015. All rights reserved, also regarding any disposal, exploitation, reproduction, editing,  
# distribution, as well as in the event of applications for industrial property rights.
# 
# Rotate mesh in stl file
# 
# 
# @author Benjamin Maier (CR/AEG)
# 
# @date 05/2016
# 

import numpy as np
from numpy import sin, cos
import random
import sys
import csv
import struct

if __name__ == '__main__':
    
    # read in arguments
    if len(sys.argv) < 6:
      print "usage: {} <fileNameIn> <fileNameOut> <x> <y> <z>\n Rotate by angles (in degrees) by x,y,z-axis by [x,y,z].".format(sys.argv[0])
      exit
    
    fileNameIn = sys.argv[1]
    fileNameOut = sys.argv[2]
    
    x = float(sys.argv[3])
    y = float(sys.argv[4])
    z = float(sys.argv[5])
    angleX = x/180.*np.pi
    angleY = y/180.*np.pi
    angleZ = z/180.*np.pi
         
          
    print "\nConvert binary STL file {} to binary stl file {}\nAngles (in degrees): [{},{},{}] (in rad): [{},{},{}]".format(fileNameIn, fileNameOut, x, y, z, angleX, angleY, angleZ)

    normalList = []
    vertexList = []

    # parse input file
    with open(fileNameIn, 'rb') as infile:
      
      header = infile.read(80)
      print "header: ",header
      
      nTriangles = struct.unpack('<i', infile.read(4))[0]
      print "nTriangles: ",nTriangles
            
      for i in range(nTriangles):
        normal = np.empty((3,1))
        normal[0] = struct.unpack('<f', infile.read(4))[0]
        normal[1] = struct.unpack('<f', infile.read(4))[0]
        normal[2] = struct.unpack('<f', infile.read(4))[0]
        vertex0 = np.empty((3,1))
        vertex1 = np.empty((3,1))
        vertex2 = np.empty((3,1))
        vertex0[0] = struct.unpack('<f', infile.read(4))[0]
        vertex0[1] = struct.unpack('<f', infile.read(4))[0]
        vertex0[2] = struct.unpack('<f', infile.read(4))[0]
        vertex1[0] = struct.unpack('<f', infile.read(4))[0]
        vertex1[1] = struct.unpack('<f', infile.read(4))[0]
        vertex1[2] = struct.unpack('<f', infile.read(4))[0]
        vertex2[0] = struct.unpack('<f', infile.read(4))[0]
        vertex2[1] = struct.unpack('<f', infile.read(4))[0]
        vertex2[2] = struct.unpack('<f', infile.read(4))[0]
        infile.read(2)
    
        normalList.append(normal)
        vertexList.append(vertex0)
        vertexList.append(vertex1)
        vertexList.append(vertex2)
      
      # compute center
      center = np.empty((3,1))
      for vertex in vertexList:
        
        center[0] = center[0]+vertex[0]
        center[1] = center[1]+vertex[1]
        center[2] = center[2]+vertex[2]
    
      print "number of vertices: ",len(vertexList)
    
      for i in range(3):
        center[i] /= len(vertexList)
            
      # transform vertices
      rotationX = np.array([[1, 0, 0], [0, np.cos(angleX), -np.sin(angleX)], [0, np.sin(angleX), np.cos(angleX)]])
      rotationY = np.array([[np.cos(angleY), 0, np.sin(angleY)], [0, 1, 0], [-np.sin(angleY), 0, np.cos(angleY)]])
      rotationZ = np.array([[np.cos(angleZ), -np.sin(angleZ), 0], [np.sin(angleZ), np.cos(angleZ), 0], [0, 0, 1]])
      
      print "rotation matrices:"
      print rotationX
      print rotationY
      print rotationZ
      
      vertexList2 = []
      
      for vertex in vertexList:
        vertex = rotationX.dot(vertex)
        vertex = rotationY.dot(vertex)
        vertex = rotationZ.dot(vertex)
      
        vertexList2.append(vertex)
        
      # compute bounding box
      minX = vertexList2[0][0]
      maxX = vertexList2[0][0]
      minY = vertexList2[0][1]
      maxY = vertexList2[0][1]
      minZ = vertexList2[0][2]
      maxZ = vertexList2[0][2]
      for vertex in vertexList2:
        minX = min(minX, vertex2[0])
        minY = min(minY, vertex2[1])
        minZ = min(minZ, vertex2[2])
        maxX = max(maxX, vertex2[0])
        maxY = max(maxY, vertex2[1])
        maxZ = max(maxZ, vertex2[2])
      
      print "bounding box: [{},{}]x[{},{}]x[{},{}]".format(minX,maxX,minY,maxY,minZ,maxZ)
      
      # write to output file  
      outfile = open (fileNameOut, "wb")
      
      # header of binary file
      outfile.write("Binary STL file")

      for i in range(80-15):
        outfile.write(struct.pack("<b", 0))
      
      # number of triangles
      outfile.write(struct.pack("<I", nTriangles))
            
      # loop over points
      for i in range(nTriangles):
        
        a = np.array([vertexList2[3*i+0][0][0], vertexList2[3*i+0][1][0], vertexList2[3*i+0][2][0]])
        b = np.array([vertexList2[3*i+1][0][0], vertexList2[3*i+1][1][0], vertexList2[3*i+1][2][0]])
        c = np.array([vertexList2[3*i+2][0][0], vertexList2[3*i+2][1][0], vertexList2[3*i+2][2][0]])
        
        # compute normal
        n = np.cross(-a+b, -a+c)
        n = n / np.linalg.norm(n)
        
        outfile.write(struct.pack("<f", n[0]))
        outfile.write(struct.pack("<f", n[1]))
        outfile.write(struct.pack("<f", n[2]))
                
        outfile.write(struct.pack("<f", a[0]))
        outfile.write(struct.pack("<f", a[1]))
        outfile.write(struct.pack("<f", a[2]))
        
        outfile.write(struct.pack("<f", b[0]))
        outfile.write(struct.pack("<f", b[1]))
        outfile.write(struct.pack("<f", b[2]))
        
        outfile.write(struct.pack("<f", c[0]))
        outfile.write(struct.pack("<f", c[1]))
        outfile.write(struct.pack("<f", c[2]))
        
        outfile.write(struct.pack("<h", 0))
        
      # close binary file
      outfile.close()
      
      print "File {} written.".format(fileNameOut)
        
      print "number of triangles: {}".format(nTriangles)

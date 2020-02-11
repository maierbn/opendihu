#!/usr/bin/env python2.7

# @file info_stl.py
# 
# Show info to stl file
# 
# 
# @author Benjamin Maier
# 
# @date 09/2017
# 

import numpy as np
from numpy import sin, cos
import random
import sys
import csv
import struct

if __name__ == '__main__':
    
    # read in arguments
    if len(sys.argv) < 2:
      print("usage: {} <fileNameIn>\n Show info.".format(sys.argv[0]))
      sys.exit(0)
    
    fileNameIn = sys.argv[1]
    
    normalList = []
    vertexList = []

    # parse input file
    with open(fileNameIn, 'rb') as infile:
      
      header = infile.read(80)
      print("Header: ",header)
      
      nTriangles = struct.unpack('<i', infile.read(4))[0]
            
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
    
      print("Number of vertices: ",len(vertexList))
      print("Number of triangles: {}".format(nTriangles))
    
      for i in range(3):
        center[i] /= len(vertexList)
            
      print("Center of gravity from vertices: {}, {}, {}".format(center[0][0], center[1][0], center[2][0]))
      
      # compute bounding box
      minX = vertexList[0][0]
      maxX = vertexList[0][0]
      minY = vertexList[0][1]
      maxY = vertexList[0][1]
      minZ = vertexList[0][2]
      maxZ = vertexList[0][2]
      for vertex in vertexList:
        minX = min(minX, vertex[0])
        minY = min(minY, vertex[1])
        minZ = min(minZ, vertex[2])
        maxX = max(maxX, vertex[0])
        maxY = max(maxY, vertex[1])
        maxZ = max(maxZ, vertex[2])
      
      print("Bounding box: [{},{}] x [{},{}] x [{},{}]".format(minX[0],maxX[0],minY[0],maxY[0],minZ[0],maxZ[0]))
      

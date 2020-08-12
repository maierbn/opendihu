#!/usr/bin/env python

# @file stl_remove_with_normals.py

import numpy as np
from numpy import sin, cos
import random
import sys
import csv
import struct

if __name__ == '__main__':
    
    # read in arguments
    if len(sys.argv) < 7:
      print "usage: {} <fileNameIn> <fileNameOut> <nx> <ny> <nz> <tol>\n Remove faces that have normal [nx,ny,nz], angle tolerance tol.".format(sys.argv[0])
      exit
    
    fileNameIn = sys.argv[1]
    fileNameOut = sys.argv[2]
    
    nx = float(sys.argv[3])
    ny = float(sys.argv[4])
    nz = float(sys.argv[5])
    tol = float(sys.argv[6])
         
          
    print "\nConvert binary STL file {} to binary stl file {}\nNormal vector: [{},{},{}]\ntolerance: {}".format(fileNameIn, fileNameOut, nx, ny, nz, tol)

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
      
      # create vertexList without faces that have normal [nx,ny,nz]
      vertexList2 = []
      normalList2 = []
      reference_normal = np.array([[nx],[ny],[nz]])
      
      print "reference normal: ",reference_normal
      
      for (i,normal) in enumerate(normalList):
        dot = normal[0]*reference_normal[0] + normal[1]*reference_normal[1] + normal[2]*reference_normal[2]
        
        if abs(np.arccos(dot)) > tol:
          normalList2.append(normal)
          vertexList2.append(vertexList[3*i+0])
          vertexList2.append(vertexList[3*i+1])
          vertexList2.append(vertexList[3*i+2])
      
      nTriangles2 = len(normalList2)
      
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
      
      print "bounding box: [{},{}]x[{},{}]x[{},{}]".format(minX,maxX,minY,maxY,minZ,maxZ)
      
      # write to output file  
      outfile = open (fileNameOut, "wb")
      
      # header of binary file
      outfile.write("Binary STL file")

      for i in range(80-15):
        outfile.write(struct.pack("<b", 0))
      
      # number of triangles
      outfile.write(struct.pack("<I", nTriangles2))
            
      # loop over points
      for i in range(nTriangles2):
        
        outfile.write(struct.pack("<f", normalList2[i][0]))
        outfile.write(struct.pack("<f", normalList2[i][1]))
        outfile.write(struct.pack("<f", normalList2[i][2]))
                
        outfile.write(struct.pack("<f", vertexList2[3*i+0][0]))
        outfile.write(struct.pack("<f", vertexList2[3*i+0][1]))
        outfile.write(struct.pack("<f", vertexList2[3*i+0][2]))
        
        outfile.write(struct.pack("<f", vertexList2[3*i+1][0]))
        outfile.write(struct.pack("<f", vertexList2[3*i+1][1]))
        outfile.write(struct.pack("<f", vertexList2[3*i+1][2]))
        
        outfile.write(struct.pack("<f", vertexList2[3*i+2][0]))
        outfile.write(struct.pack("<f", vertexList2[3*i+2][1]))
        outfile.write(struct.pack("<f", vertexList2[3*i+2][2]))
        
        outfile.write(struct.pack("<h", 0))
        
      # close binary file
      outfile.close()
      
      print "File {} written.".format(fileNameOut)
        
      print "Number of removed triangles: {} of {}, remaining: {}".format(nTriangles-nTriangles2, nTriangles, nTriangles2)

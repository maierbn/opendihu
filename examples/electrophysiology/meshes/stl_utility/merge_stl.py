#!/usr/bin/env python

# @file center_stl.py
# 
# (c) Robert Bosch GmbH 2015. All rights reserved, also regarding any disposal, exploitation, reproduction, editing,  
# distribution, as well as in the event of applications for industrial property rights.
# 
# Translate mesh in stl file according to center of mass of vertices.
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
    if len(sys.argv) < 3:
      print "usage: {} <fileNameIn1> <fileNameIn2> <fileNameOut>\n Merge files.".format(sys.argv[0])
      exit
    
    fileNameIn1 = sys.argv[1]
    fileNameIn2 = sys.argv[2]
    fileNameOut = sys.argv[3]
    
    print "\Merge binary STL files {}, {} to binary stl file {}".format(fileNameIn1, fileNameIn2, fileNameOut)

    normalList = []
    vertexList = []

    # parse input file
    with open(fileNameIn1, 'rb') as infile1:
      with open(fileNameIn2, 'rb') as infile2:
        
        header = infile1.read(80)
        header = infile2.read(80)
        print "header: ",header
        
        nTriangles1 = struct.unpack('<i', infile1.read(4))[0]
        nTriangles2 = struct.unpack('<i', infile2.read(4))[0]
        nTriangles = nTriangles1+nTriangles2
        print "nTriangles: {} + {} = {}".format(nTriangles1,nTriangles2,nTriangles)
              
        # read file 1
        for i in range(nTriangles1):
          normal = np.empty((3,1))
          normal[0] = struct.unpack('<f', infile1.read(4))[0]
          normal[1] = struct.unpack('<f', infile1.read(4))[0]
          normal[2] = struct.unpack('<f', infile1.read(4))[0]
          vertex0 = np.empty((3,1))
          vertex1 = np.empty((3,1))
          vertex2 = np.empty((3,1))
          vertex0[0] = struct.unpack('<f', infile1.read(4))[0]
          vertex0[1] = struct.unpack('<f', infile1.read(4))[0]
          vertex0[2] = struct.unpack('<f', infile1.read(4))[0]
          vertex1[0] = struct.unpack('<f', infile1.read(4))[0]
          vertex1[1] = struct.unpack('<f', infile1.read(4))[0]
          vertex1[2] = struct.unpack('<f', infile1.read(4))[0]
          vertex2[0] = struct.unpack('<f', infile1.read(4))[0]
          vertex2[1] = struct.unpack('<f', infile1.read(4))[0]
          vertex2[2] = struct.unpack('<f', infile1.read(4))[0]
          infile1.read(2)
      
          normalList.append(normal)
          vertexList.append(vertex0)
          vertexList.append(vertex1)
          vertexList.append(vertex2)
        
        # read file 2
        for i in range(nTriangles2):
          normal = np.empty((3,1))
          normal[0] = struct.unpack('<f', infile2.read(4))[0]
          normal[1] = struct.unpack('<f', infile2.read(4))[0]
          normal[2] = struct.unpack('<f', infile2.read(4))[0]
          vertex0 = np.empty((3,1))
          vertex1 = np.empty((3,1))
          vertex2 = np.empty((3,1))
          vertex0[0] = struct.unpack('<f', infile2.read(4))[0]
          vertex0[1] = struct.unpack('<f', infile2.read(4))[0]
          vertex0[2] = struct.unpack('<f', infile2.read(4))[0]
          vertex1[0] = struct.unpack('<f', infile2.read(4))[0]
          vertex1[1] = struct.unpack('<f', infile2.read(4))[0]
          vertex1[2] = struct.unpack('<f', infile2.read(4))[0]
          vertex2[0] = struct.unpack('<f', infile2.read(4))[0]
          vertex2[1] = struct.unpack('<f', infile2.read(4))[0]
          vertex2[2] = struct.unpack('<f', infile2.read(4))[0]
          infile2.read(2)
      
          normalList.append(normal)
          vertexList.append(vertex0)
          vertexList.append(vertex1)
          vertexList.append(vertex2)
        
        print "number of vertices: ",len(vertexList)
      
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
        outfile.write(struct.pack("<I", nTriangles))
              
        # loop over points
        for i in range(nTriangles):
          
          outfile.write(struct.pack("<f", normalList[i][0]))
          outfile.write(struct.pack("<f", normalList[i][1]))
          outfile.write(struct.pack("<f", normalList[i][2]))
                  
          outfile.write(struct.pack("<f", vertexList[3*i+0][0]))
          outfile.write(struct.pack("<f", vertexList[3*i+0][1]))
          outfile.write(struct.pack("<f", vertexList[3*i+0][2]))
          
          outfile.write(struct.pack("<f", vertexList[3*i+1][0]))
          outfile.write(struct.pack("<f", vertexList[3*i+1][1]))
          outfile.write(struct.pack("<f", vertexList[3*i+1][2]))
          
          outfile.write(struct.pack("<f", vertexList[3*i+2][0]))
          outfile.write(struct.pack("<f", vertexList[3*i+2][1]))
          outfile.write(struct.pack("<f", vertexList[3*i+2][2]))
          
          outfile.write(struct.pack("<h", 0))
          
        # close binary file
        outfile.close()
        
        print "File {} written.".format(fileNameOut)
        
      print "Total number of triangles: {}".format(nTriangles)

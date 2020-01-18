#!/usr/bin/env python

# @file rotations.py
# 
# (c) Robert Bosch GmbH 2015. All rights reserved, also regarding any disposal, exploitation, reproduction, editing,  
# distribution, as well as in the event of applications for industrial property rights.
# 
# Convert .msh files to binary .stl files
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
    if len(sys.argv) < 4:
      print "usage: {} <fileNameIn> <fileNameOutBinary> <fileNameOutAscii>".format(sys.argv[0])
      exit
    
    fileNameIn = sys.argv[1]
    fileNameOut = sys.argv[2]
    fileNameOutAscii = sys.argv[3]
    
    print "\nConvert file {}\n * to binary stl: {}\n * to ascii stl: {}\n".format(fileNameIn, fileNameOut, fileNameOutAscii)

    # parse file
    with open(fileNameIn) as csvfile:
      reader = csv.reader(csvfile, delimiter=' ')
      
      vertexCount = 0
      edgeCount = 0
      faceCount = 0
      
      vertexList = []
      edgeList = []
      faceList = []
      
      i = 0
      
      currentMode = ''
      
      # read in lists
      for row in reader:
        if len(row) >= 2:
          if row[0] == 'Vertex_Count:':
            vertexCount = int(row[1])
            vertexList = np.empty((vertexCount, 4))
            currentMode = 'vertex'
            i = 0
            
          elif row[0] == 'Edge_Count:':
            edgeCount = int(row[1])
            edgeList = np.empty((edgeCount, 2))
            currentMode = 'edge'
            i = 0
            
          elif row[0] == 'Face_Count:':
            faceCount = int(row[1])
            faceList = np.empty((faceCount, 3))
            currentMode = 'face'
            i = 0
          
          elif currentMode == 'vertex':
            for j in range(4):
              vertexList[i][j] = float(row[j])
            i += 1
              
          elif currentMode == 'edge':
            for j in range(2):
              edgeList[i][j] = int(row[j])
            i += 1
              
          elif currentMode == 'face':
            for j in range(3):
              faceList[i][j] = int(row[j])
            i += 1
            
      print "File {} loaded".format(fileNameIn)
            
      # write to output file  
      outfile = open (fileNameOut, "wb")
      
      # open ascii file
      outfileAscii = open (fileNameOutAscii, "w")
      outfileAscii.write("solid mesh\n")
      
      # header of binary file
      #bytesList = 80*[0]
      #print bytesList
      #outfile.write(bytearray(bytesList))
      outfile.write("Binary STL file")

      for i in range(80-15):
        outfile.write(struct.pack("<b", 0))
      
      # number of triangles
      outfile.write(struct.pack("<I", faceCount))
            
            
      minPoint = np.ones((3))*10000.
      maxPoint = np.zeros((3))
      
      # loop over points
      for i in range(faceCount):
        
        # get points
        a = vertexList[faceList[i][0]][0:3]
        b = vertexList[faceList[i][1]][0:3]
        c = vertexList[faceList[i][2]][0:3]
        
        minPoint[0] = min(a[0], minPoint[0])
        minPoint[1] = min(a[1], minPoint[1])
        minPoint[2] = min(a[2], minPoint[2])
        
        maxPoint[0] = max(a[0], maxPoint[0])
        maxPoint[1] = max(a[1], maxPoint[1])
        maxPoint[2] = max(a[2], maxPoint[2])
        
        minPoint[0] = min(b[0], minPoint[0])
        minPoint[1] = min(b[1], minPoint[1])
        minPoint[2] = min(b[2], minPoint[2])
        
        maxPoint[0] = max(b[0], maxPoint[0])
        maxPoint[1] = max(b[1], maxPoint[1])
        maxPoint[2] = max(b[2], maxPoint[2])
        
        minPoint[0] = min(c[0], minPoint[0])
        minPoint[1] = min(c[1], minPoint[1])
        minPoint[2] = min(c[2], minPoint[2])
        
        maxPoint[0] = max(c[0], maxPoint[0])
        maxPoint[1] = max(c[1], maxPoint[1])
        maxPoint[2] = max(c[2], maxPoint[2])
        
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
        
        # write to ascii file
        outfileAscii.write("facet normal {} {} {}\n".format(n[0], n[1], n[2]))
        
        outfileAscii.write("\touter loop\n")
        outfileAscii.write("\t\tvertex {} {} {}\n".format(a[0], a[1], a[2]))
        outfileAscii.write("\t\tvertex {} {} {}\n".format(b[0], b[1], b[2]))
        outfileAscii.write("\t\tvertex {} {} {}\n".format(c[0], c[1], c[2]))
        outfileAscii.write("\tendloop\n")
        outfileAscii.write("endfacet\n")
        
      # close binary file
      outfile.close()
      
      # close ascii file
      outfileAscii.write("endsolid mesh\n")
      outfileAscii.close()
      
        
      print "Files {}, {} written.".format(fileNameOut, fileNameOutAscii)
        
      print "Statistics: vertexCount: {}, edgeCount: {}, faceCount: {}".format(vertexCount, edgeCount, faceCount)
      #print "vertexList: {}\n, edgeList: {}\n, faceList: {}".format(vertexList, edgeList, faceList)
      print "Bounding box: min: {}, max: {}".format(minPoint, maxPoint)

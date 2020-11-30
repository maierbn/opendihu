#!/usr/bin/env ../../../../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# The functions in this script will be called from the C++ code from parallel_fiber_estimation. They write some debugging information to files.

import sys, os
import numpy as np
import stl
from stl import mesh
import vtk

import stl_create_rings

def output_points(filename, rankNo, level, points, size):
  
  with_vtk = False
  triangles = []
  #print("> output_points(filename={}, rankNo={}, level={}, n points: {}, size={})".format(filename, rankNo, level, len(points), size))
  
  if with_vtk:
    try:
      # setup points and vertices
      vtk_points = vtk.vtkPoints()
    except:
      pass
  
  factor = 1.0
  for p in points:
    point = np.array([p[0], p[1], p[2]])
    
    # add point to vtk data set
    if with_vtk:
      vtk_points.InsertNextPoint(p[0],p[1],p[2])
  
    # add triangle to stl dataset
    stl_create_rings.create_point_marker(point, triangles, size*factor)
    
  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  if level != -1:
    if not os.path.exists("out/level_{}".format(level)):
      os.makedirs("out/level_{}".format(level),0o755)
    outfile = "out/level_{}/{}.{}.{}.stl".format(level, filename[0:2], rankNo, filename[2:])
  else:
    outfile = "{}.stl".format(filename)
    
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))
  
  if with_vtk:
    try:
      polydata = vtk.vtkPolyData()
      polydata.SetPoints(vtk_points)
      polydata.Modified()

      writer = vtk.vtkXMLPolyDataWriter()
      writer.SetFileName(outfile+".vtp")
      writer.SetInputData(polydata)
      writer.Write()
    except:
      print("writing vtp file {} failed".format(filename))
    
  print("> output_points(filename={}, rankNo={}, level={}, n points: {}, size={}) done".format(filename, rankNo, level, len(points), size))

def output_streamline(filename, rankNo, level, points, size):
  
  triangles = []

  print("> output_streamline(filename={}, rankNo={}, level={}, n points: {}, size={})".format(filename, rankNo, level, len(points), size))
  
  try:  
     # setup points and vertices
    vtk_points = vtk.vtkPoints()
    vtk_lines = vtk.vtkCellArray()
  except:
    pass
  
  factor = 1.0
  line_no = 0

  previous_point = None
  for p in points:
    point = np.array([p[0], p[1], p[2]])
    if np.linalg.norm(point) < 1e-3:
      continue
    if previous_point is not None:
      triangles.append([previous_point, point, 0.5*(previous_point+point)])      
            
      try:
        vtk_points.InsertNextPoint(previous_point[0], previous_point[1], previous_point[2])
        vtk_points.InsertNextPoint(point[0], point[1], point[2])
        
        vtk_line = vtk.vtkLine()
        vtk_line.GetPointIds().SetId(0, 2*line_no + 0)
        vtk_line.GetPointIds().SetId(1, 2*line_no + 1)
        vtk_lines.InsertNextCell(vtk_line)
        line_no += 1
      except:
        print("Error in creating vtk dataset in output_streamline({})".format(filename))
      
    previous_point = point
    
    stl_create_rings.create_point_marker(point, triangles, size*factor)
    #factor *= 1.1
    #if factor > 3:
    #  factor = 3.0

  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  if not os.path.exists("out/level_{}".format(level)):
    os.makedirs("out/level_{}".format(level),0o755)
  outfile = "out/level_{}/{}.{}.{}.stl".format(level, filename[0:2], rankNo, filename[2:])
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))
  
  try:
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetLines(vtk_lines)

    polydata.Modified()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outfile+".vtp")
    writer.SetInputData(polydata)
    writer.Write()
  except:
    print("writing vtp file {} failed".format(filename))
    
  print("> output_streamline(filename={}, rankNo={}, level={}, n points: {}, size={}) done".format(filename, rankNo, level, len(points), size))

def output_streamlines(filename, rankNo, level, streamlines, size):
  
  triangles = []
  
  try:
     # setup points and vertices
    vtk_points = vtk.vtkPoints()
    vtk_lines = vtk.vtkCellArray()
  except:
    pass
  
  factor = 1.0
  line_no = 0
  
  for points in streamlines:
    previous_point = None
    
    #print("output_streamlines, streamline: {}".format(points))
    for p in points:
      point = np.array([p[0], p[1], p[2]])
      if np.linalg.norm(point) < 1e-3:
        continue
      if previous_point is not None:
        triangles.append([previous_point, point, 0.5*(previous_point+point)])
          
        try:
          vtk_points.InsertNextPoint(previous_point[0], previous_point[1], previous_point[2])
          vtk_points.InsertNextPoint(point[0], point[1], point[2])
          
          vtk_line = vtk.vtkLine()
          vtk_line.GetPointIds().SetId(0, 2*line_no + 0)
          vtk_line.GetPointIds().SetId(1, 2*line_no + 1)
          vtk_lines.InsertNextCell(vtk_line)
        except:
          print("Error in creating vtk dataset in output_streamlines({})".format(filename))
          
        line_no += 1
      previous_point = point
      
      stl_create_rings.create_point_marker(point, triangles, size*factor)
      #factor *= 1.1
      #if factor > 3:
      #  factor = 3.0

  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  if not os.path.exists("out/level_{}".format(level)):
    os.makedirs("out/level_{}".format(level),0o755)
  outfile = "out/level_{}/{}.{}.{}.stl".format(level, filename[0:2], rankNo, filename[2:])
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))
  
  try:
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetLines(vtk_lines)

    polydata.Modified()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outfile+".vtp")
    writer.SetInputData(polydata)
    writer.Write()
  except:
    print("writing vtp file {} failed".format(filename))

def output_rings(filename, rankNo, level, rings, size):
  
  triangles = []
  print("> output_rings(filename={}, rankNo={}, level={}, n rings: {}, size={})".format(filename, rankNo, level, len(rings), size))
  
  # setup points and vertices
  try:
    vtk_points = vtk.vtkPoints()
    vtk_lines = vtk.vtkCellArray()
  except Exception as error:
    print("Error in setup for output_rings {}: {}".format(filename, error))
  
  
  factor = 1.0
  line_no = 0
  
  for points in rings:
    previous_point = None
    first_point = None
    
    for p in points:
      point = np.array([p[0], p[1], p[2]])
      if np.linalg.norm(point) < 1e-3:
        continue
      if previous_point is None:
        first_point = point
      else:
        triangles.append([previous_point, point, 0.5*(previous_point+point)])
        
        try:
          vtk_points.InsertNextPoint(previous_point[0], previous_point[1], previous_point[2])
          vtk_points.InsertNextPoint(point[0], point[1], point[2])
          
          vtk_line = vtk.vtkLine()
          vtk_line.GetPointIds().SetId(0, 2*line_no+0)
          vtk_line.GetPointIds().SetId(1, 2*line_no+1)
          vtk_lines.InsertNextCell(vtk_line)
        except Exception as error:
          print("Error in creating vtk dataset in output_rings({}): {}".format(filename, error))
          
        line_no += 1
      previous_point = point
      
      stl_create_rings.create_point_marker(point, triangles, size*factor)
      #factor *= 1.1
      #if factor > 3:
      #  factor = 3.0

    # close loop (not for boundary points on faces)
    if previous_point is not None:
      triangles.append([previous_point, first_point, 0.5*(previous_point+first_point)])
    
      try:
        vtk_points.InsertNextPoint(previous_point[0], previous_point[1], previous_point[2])
        vtk_points.InsertNextPoint(first_point[0], first_point[1], first_point[2])
        
        vtk_line = vtk.vtkLine()
        vtk_line.GetPointIds().SetId(0, 2*line_no+0)
        vtk_line.GetPointIds().SetId(1, 2*line_no+1)
        vtk_lines.InsertNextCell(vtk_line)
      except Exception as error:
        print("Error in creating vtk dataset in output_rings({}): {}".format(filename, error))
        

  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  print("output_rings({}, {}, {}, ...)".format(filename, rankNo, level))

  if not os.path.exists("out/level_{}".format(level)):
    print("path does not exist")
    os.makedirs("out/level_{}".format(level),0o755)
    print("path created")
  outfile = "out/level_{}/{}.{}.{}.stl".format(level, filename[0:2], rankNo, filename[2:])
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))
    
  try:
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetLines(vtk_lines)

    polydata.Modified()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outfile+".vtp")
    writer.SetInputData(polydata)
    writer.Write()
  except Exception as error:
    print("writing vtp file {} failed: {}".format(filename, error))
    
  print("> output_rings(filename={}, rankNo={}, n rings: {}, size={}) done".format(filename, rankNo, level, len(rings), size))

def output_boundary_points(filename, rankNo, level, points, size):
  
  triangles = []
  
  print("> output_boundary_points(filename={}, rankNo={}, level={}, n points: {}, size={})".format(filename, rankNo, level, len(points), size))
  
  try:
    # setup points and vertices
    vtk_points = vtk.vtkPoints()
    vtk_lines = vtk.vtkCellArray()
  except Exception as error:
    print("Error in setup for output_boundary_points, {}: {}".format(filename, error))
  
  
  # data structure:
  # std::array<std::vector<std::vector<Vec3>>,4>
  # list [<face0>, <face1>, <face2>, <face3>]
  # <face> = [<level0>, <level1>, ...]
  # <level> = [<point0>, <point1>, ...]
  
  factor = 1.0
  line_no = 0
  for face_points in points:
    for zlevel_points in face_points:
      previous_point = None
      first_point = None
      for p in zlevel_points:
        point = np.array([p[0], p[1], p[2]])
        if np.linalg.norm(point) < 1e-3:
          continue
        if previous_point is None:
          first_point = point
        else:
          triangles.append([previous_point, point, 0.5*(previous_point+point)])
            
          try:
            id1 = vtk_points.InsertNextPoint(previous_point[0], previous_point[1], previous_point[2])
            id2 = vtk_points.InsertNextPoint(point[0], point[1], point[2])
            
            vtk_line = vtk.vtkLine()
            vtk_line.GetPointIds().SetId(0, id1)
            vtk_line.GetPointIds().SetId(1, id2)
            vtk_lines.InsertNextCell(vtk_line)
            line_no += 1
          except Exception as error:
            print("Error in creating vtk dataset in output_boundary_points({}): {}".format(filename, error))
            
        previous_point = point
        
        stl_create_rings.create_point_marker(point, triangles, size*factor)
        #factor *= 1.1
        #if factor > 3:
        #  factor = 3.0        
      
      # close loop (not for boundary points on faces)
      #if previous_point is not None:
      #  triangles.append([previous_point, first_point, 0.5*(previous_point+first_point)])

  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
    #for j in range(3):
      #print("set (",i,",",j,")=",f[j]," (=",stl_mesh.vectors[i][j],")"
      
  #out_mesh.update_normals()

  if not os.path.exists("out/level_{}".format(level)):
    os.makedirs("out/level_{}".format(level),0o755)
  outfile = "out/level_{}/{}.{}.{}.stl".format(level, filename[0:2], rankNo, filename[2:])
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))
  
  try:
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetLines(vtk_lines)

    polydata.Modified()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outfile+".vtp")
    writer.SetInputData(polydata)
    writer.Write()
  except Exception as error:
    print("writing vtp file {} failed: {}".format(filename, error))
    
  print("> output_boundary_points(filename={}, rankNo={}, level={}, n points: {}, size={}) done".format(filename, rankNo, level, len(points), size))

def output_triangles(filename, triangles):
  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  outfile = "out/{}.stl".format(filename)
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))

def output_ghost_elements(filename, rankNo, level, point_values, n_elements, size):
  
  triangles = []
  
  try:
    # setup points and vertices
    vtk_points = vtk.vtkPoints()
    vtk_boxes = vtk.vtkCellArray()
  except:
    pass
  
  n_points = (int)(len(point_values)/3)
  
  n_nodes = [n_elements[0]+1, n_elements[1]+1, n_elements[2]+1]
  
  #print("output_ghost_elements, filename {}, n_elements: {}, n points: {}, n_nodes: {}, n_points: {}".format(filename, n_elements, len(point_values), n_nodes, n_points))
  #print("point_values: {}".format(point_values))
  #print("n_points: {}".format(n_points))
  #print("n_elements: {}".format(n_elements))
  
  factor = 1.0
  box_no = 0
  for z in range(n_elements[2]):
    for y in range(n_elements[1]):
      for x in range(n_elements[0]):
        p = list()
        
        for k in range(2):
          for j in range(2):
            for i in range(2):
              index = (z+k) * n_nodes[0]*n_nodes[1] + (y+j) * n_nodes[0] + (x+i)
              #print("index: {}".format(index))
              p0 = np.array([point_values[index], point_values[n_points+index], point_values[2*n_points+index]])
              p.append(p0)
        #print("x: {}, y: {}, z: {}, points: {}".format(x,y,z,p))
              
        triangles += [
          [p[0],p[3],p[1]],[p[0],p[2],p[3]],  # bottom
          [p[4],p[5],p[7]],[p[4],p[7],p[6]],  # top
          [p[0],p[1],p[5]],[p[0],p[5],p[4]],  # front
          [p[2],p[7],p[3]],[p[2],p[6],p[7]],  # back
          [p[2],p[0],p[4]],[p[2],p[4],p[6]],  # left
          [p[1],p[3],p[7]],[p[1],p[7],p[5]]  # right
        ]
        
        try:
          for i in [0,1,3,2,4,5,7,6]:
            point = p[i]
            vtk_points.InsertNextPoint(point[0], point[1], point[2])
          
          vtk_box = vtk.vtkLine()
          for i in range(8):
            vtk_box.GetPointIds().SetId(i, 8*box_no + i)
          vtk_boxes.InsertNextCell(8)
          
        except:
          print("Error in creating vtk dataset in output_ghost_elements({})".format(filename))
        box_no += 1
  #---------------------------------------
  # Create the mesh
  out_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
  for i, f in enumerate(triangles):
    out_mesh.vectors[i] = f
  #out_mesh.update_normals()

  if not os.path.exists("out/level_{}".format(level)):
    os.makedirs("out/level_{}".format(level),0o755)
  outfile = "out/level_{}/{}.{}.{}.stl".format(level, filename[0:2], rankNo, filename[2:])
  #out_mesh.save(outfile, mode=stl.Mode.ASCII)
  out_mesh.save(outfile)
  print("saved {} triangles to \"{}\"".format(len(triangles),outfile))

  try:
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetPolys(vtk_boxes)

    polydata.Modified()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outfile+".vtp")
    writer.SetInputData(polydata)
    writer.Write()
  except:
    print("writing vtp file {} failed".format(filename))



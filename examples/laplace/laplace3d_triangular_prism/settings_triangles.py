import numpy as np
import sys


import plot_quadrangulation_schemes

# 3D laplace problem

nx = 2
ny = 2
nz = 1

nx = 3
ny = 3
nz = 8

quadratic = True
mx = nx+1
my = ny+1
mz = nz+1
if quadratic:
  mx = 2*nx+1
  my = 2*ny+1
  mz = 2*nz+1

# node positions
node_positions = []

for k in range(mz):
  z = k
  if quadratic:
    z = k/2
  
  if quadratic and False:
    node_positions += [
      [-0.5,-0.5,z], [-0.25,-0.75,z], [0,-1,z], [0.25,-0.75,z], [0.5,-0.5,z],
      [-0.75,-0.25,z], [-0.25,-0.5,z], [0,-0.5,z], [0.25,-0.5,z], [0.75,-0.25,z],
      [-1,0,z], [-0.5,0,z], [0,0,z], [0.5,0.0,z], [1,0,z],
      [-0.75,0.25,z], [-0.25,0.5,z], [0,0.5,z], [0.25,0.5,z], [0.75,0.25,z],
      [-0.5,0.5,z], [-0.25,0.75,z], [0,1,z], [0.25,0.75,z], [0.5,0.5,z]]
  else:  
    n_grid_points_x = mx
    n_grid_points_y = my
    parametric_space_shape = 3
    node_positions_z_plane = plot_quadrangulation_schemes.create_reference_domain_quadrangulation(n_grid_points_x, n_grid_points_y, parametric_space_shape)

    r = np.sqrt(2)
    alpha = np.pi/4
    r = 3
    node_positions += [[r*(np.sin(alpha)*p[0] + np.cos(alpha)*p[1]), r*(-np.cos(alpha)*p[0] + np.sin(alpha)*p[1]), z] for p in node_positions_z_plane]

# Dirichlet boundary conditions
bc = {}

for j in range(my):
  for i in range(mx):
    # x = 0 plane
    x = i/(mx-1)
    y = j/(my-1)
    index = j*mx + i
  
    # z- plane (bottom)
    bc[index] = np.sin(x*np.pi)
    bc[index] = 1.0
    
    # z+ plane (top)
    index += (mz-1)*mx*my
    bc[index] = np.sin(y*np.pi) + 2.0
    bc[index] = 2.0

rank_no = (int)(sys.argv[-2])

config = {
  "scenarioName": "",
  "mappingsBetweenMeshesLogFile": None,
  "logFormat":                      "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "FiniteElementMethod" : {
    "nElements":  [nx,ny,nz],
    "nRanks":     [1,1,1],
    "hasTriangleCorners": True,
    "inputMeshIsGlobal": True,
    "nodePositions": node_positions,
    "outputInterval": 1,
    "slotName": None,
    
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions":   [],
    "prefactor": 1,
    
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 10000,
    "dumpFormat": "matlab",
    "dumpFilename": "out/a",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},      
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"}
    ]
  },
}

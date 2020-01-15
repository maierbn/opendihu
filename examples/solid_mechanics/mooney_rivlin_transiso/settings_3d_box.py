# Mooney-Rivlin transversely isotropic with parameters from Heidlauf 2013

import numpy as np
import sys, os

# number of elements
nx = 2    # 2
ny = 1    # 2
nz = 5    # 5

# boundary conditions (for quadratic elements)
dirichlet_bc = {}

xpos = 0.0
ypos = 0.0
zpos = 0.0

own_rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# bottom plane
for j in range(0,2*ny+1):
  for i in range(0,2*nx+1):
    dirichlet_bc[j*(2*nx+1) + i] = [None,None,zpos]

if False:
  # left plane
  for k in range(0,2*nz+1):
    for j in range(0,2*ny+1):
      dirichlet_bc[k*(2*nx+1)*(2*ny+1) + j*(2*nx+1)] = [xpos,None,None]

  # front plane
  for k in range(0,2*nz+1):
    for i in range(0,2*nx+1):
      dirichlet_bc[k*(2*nx+1)*(2*ny+1) + i] = [None,ypos,None]

  # vertical edge
  for k in range(0,2*nz+1):
    dirichlet_bc[k*(2*nx+1)*(2*ny+1)] = [xpos,ypos,None]

if True:
  # horizontal edge
  for i in range(0,2*nx+1):
    dirichlet_bc[i] = [None,ypos,zpos]

if False:
  # horizontal edge
  for j in range(0,2*ny+1):
    dirichlet_bc[j*(2*nx+1)] = [xpos,None,zpos]

# corner
dirichlet_bc[0] = [xpos,ypos,zpos]
dirichlet_bc[1] = [None,ypos,zpos]

neumann_bc = [{"element": (nz-1)*nx*ny + j*nx + i, "constantVector": [1e-1,0,0], "face": "2+"} for j in range(ny) for i in range(nx)]

#dirichlet_bc = {}
#neumann_bc = []

# fiber directions
fiber_meshes = {}
fiber_mesh_names = []

n_elements = 0
if own_rank_no==0:
  n_elements =  2*int(np.ceil(nz/2.))+1    # 7
elif own_rank_no == 1:
  n_elements =  2*(nz-int(np.ceil(nz/2.)))  # 4

for j in range(ny*2+2):
  for i in range(nx):
    fiber_no = j*nx + i
    
    x = 0.5 + i
    y = -2.0 + 0.2 + j
    angle = 20./180.*np.pi
    
    node_positions = []
    for z in range(2*nz+1):
      if own_rank_no==0 and z > 2*int(np.ceil(nz/2.)):
        continue
      
      if own_rank_no==1 and z < 2*int(np.ceil(nz/2.)):
        continue
      
      h = 0.5*z
      x_pos = x
      y_pos = y + np.sin(angle)*h
      z_pos = 0.2 + np.cos(angle)*h
      node_positions.append([x_pos,y_pos,z_pos])
    
    mesh_name = "fiber{}".format(fiber_no)
    fiber_mesh_names.append(mesh_name)
    
    fiber_meshes[mesh_name] = {
      "nodePositions": node_positions,
      "nElements": [n_elements],
      "inputMeshIsGlobal": False,
      "nRanks": [2],
    }
    
config = {
  "scenarioName": "3d_box",
  "Meshes": fiber_meshes,
  "HyperelasticitySolver": {
    "durationLogKey": "nonlinear",
    
    #"materialParameters": [6.352e-10, 3.627e-10, 2.756e-5, 43.373],  # c1, c2, b1, d1
    "materialParameters": [2.0, 3, 4, 5],  # c1, c2, b1, d1
    "displacementsScalingFactor": 1,
    
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": False,
    "useNumericJacobian": True,   # Only works in parallel execution. If both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables": False,   # extra output of matlab vectors, x,r, jacobian matrix
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    "fiberMeshNames": fiber_mesh_names,   # fiber meshes that will be used to determine the fiber direction
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": [nx, ny, nz],
    
    # solver
    "relativeTolerance": 1e-10,
    "solverType": "preonly",          # cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "lu",
    "maxIterations": 1e4,
    "dumpFilename": "out/m",
    "dumpFormat": "matlab",   # default, ascii, matlab
    
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    
    "OutputWriter" : [   # output files for displacements function space (quadratic elements)
      {"format": "Paraview", "outputInterval": 1, "filename": "out/u", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out/u", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ],
    "pressure": {   # output files for pressure function space (linear elements)
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
        {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
      ]
    }
  },
}

print(config)

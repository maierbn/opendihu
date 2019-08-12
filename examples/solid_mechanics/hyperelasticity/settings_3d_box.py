import numpy as np
import sys, os

# number of elements
nx = 4
ny = 3
nz = 5

# boundary conditions (for quadratic elements)
dirichlet_bc = {}

xpos = 1.0
ypos = 2.0
zpos = 3.0

# left plane
for k in range(0,2*nz+1):
  for j in range(0,2*ny+1):
    dirichlet_bc[k*(2*nx+1)*(2*ny+1) + j*(2*nx+1)] = [xpos,np.nan,np.nan]

# front plane
for k in range(0,2*nz+1):
  for i in range(0,2*nx+1):
    dirichlet_bc[k*(2*nx+1)*(2*ny+1) + i] = [np.nan,ypos,np.nan]

# bottom plane
for j in range(0,2*ny+1):
  for i in range(0,2*nx+1):
    dirichlet_bc[j*(nx+1) + i] = [np.nan,np.nan,zpos]

# vertical edge
for k in range(0,2*nz+1):
  dirichlet_bc[k*(nx+1)*(ny+1)] = [xpos,ypos,np.nan]

# horizontal edge
for i in range(0,2*nx+1):
  dirichlet_bc[i] = [np.nan,ypos,zpos]

# horizontal edge
for j in range(0,2*ny+1):
  dirichlet_bc[j*(2*nx+1)] = [xpos,np.nan,zpos]

# corner
dirichlet_bc[0] = [xpos,ypos,zpos]

neumann_bc = [{"element": k*nx*ny + j*nx + nx-1, "constantVector": [+0.1,+0.6,3.0], "face": "0+"} for k in range(nz) for j in range(ny)]

#dirichlet_bc = {}
#neumann_bc = []

config = {
  "scenarioName": "3d_box",
  "QuasiStaticHyperelasticitySolver": {
    "c0": 1.0,       # dummy value
    "c1": 1.0,    # dummy value
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": True,
    "useNumericJacobian": False,   # only works with non-nested matrices
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": [nx, ny, nz],
    
    # solver
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "none",
    "maxIterations": 1e4,
    
    "dumpFilename": "out/m",
    "dumpFormat": "ascii",   # default, ascii, matlab
    
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/out", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out/out", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}

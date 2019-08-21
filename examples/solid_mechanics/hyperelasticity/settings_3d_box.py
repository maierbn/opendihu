import numpy as np
import sys, os

# number of elements
nx = 1
ny = 1
nz = 2

# boundary conditions (for quadratic elements)
dirichlet_bc = {}

xpos = 0.0
ypos = 0.0
zpos = 0.0

# bottom plane
for j in range(0,2*ny+1):
  for i in range(0,2*nx+1):
    dirichlet_bc[j*(2*nx+1) + i] = [np.nan,np.nan,zpos]

if False:
  # left plane
  for k in range(0,2*nz+1):
    for j in range(0,2*ny+1):
      dirichlet_bc[k*(2*nx+1)*(2*ny+1) + j*(2*nx+1)] = [xpos,np.nan,np.nan]

  # front plane
  for k in range(0,2*nz+1):
    for i in range(0,2*nx+1):
      dirichlet_bc[k*(2*nx+1)*(2*ny+1) + i] = [np.nan,ypos,np.nan]

  # vertical edge
  for k in range(0,2*nz+1):
    dirichlet_bc[k*(2*nx+1)*(2*ny+1)] = [xpos,ypos,np.nan]

  # horizontal edge
  for i in range(0,2*nx+1):
    dirichlet_bc[i] = [np.nan,ypos,zpos]

  # horizontal edge
  for j in range(0,2*ny+1):
    dirichlet_bc[j*(2*nx+1)] = [xpos,np.nan,zpos]

# corner
dirichlet_bc[0] = [xpos,ypos,zpos]
dirichlet_bc[1] = [np.nan,ypos,zpos]

neumann_bc = [{"element": k*nx*ny + j*nx + nx-1, "constantVector": [0,0,1], "face": "2+"} for k in range(nz) for j in range(ny)]

#dirichlet_bc = {}
#neumann_bc = []

config = {
  "scenarioName": "3d_box",
  "QuasiStaticHyperelasticitySolver": {
    "c1": 1.0,       # dummy value
    "c2": 0.0,    # dummy value
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": False,
    "useNumericJacobian": True,   # only works with non-nested matrices
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": [nx, ny, nz],
    
    # solver
    "relativeTolerance": 1e-15,
    "solverType": "gmres",          # cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "none",
    "maxIterations": 1e4,
    
    "dumpFilename": "out/m",
    "dumpFormat": "matlab",   # default, ascii, matlab
    
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/out", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out/out", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}

import numpy as np
import sys, os

nx = 5
ny = 4
nz = 3

# boundary conditions (for linear elements)
dirichlet_bc = {}

# left plane
for k in range(0,nz+1):
  for j in range(0,ny+1):
    dirichlet_bc[k*(nx+1)*(ny+1) + j*(nx+1)] = [0.0,np.nan,np.nan]

# front plane
for k in range(0,nz+1):
  for i in range(0,nx+1):
    dirichlet_bc[k*(nx+1)*(ny+1) + i] = [np.nan,0.0,np.nan]

# bottom plane
for j in range(0,ny+1):
  for i in range(0,nx+1):
    dirichlet_bc[j*(nx+1) + i] = [np.nan,np.nan,0.0]

# vertical edge
for k in range(0,nz+1):
  dirichlet_bc[k*(nx+1)*(ny+1)] = [0.0,0.0,np.nan]

# horizontal edge
for i in range(0,nx+1):
  dirichlet_bc[i] = [np.nan,0.0,0.0]

# horizontal edge
for j in range(0,ny+1):
  dirichlet_bc[j*(nx+1)] = [0.0,np.nan,0.0]

# corner
dirichlet_bc[0] = [0.0,0.0,0.0]

neumann_bc = [{"element": k*nx*ny + j*nx + nx-1, "constantVector": [+0.1,+0.6,3.0], "face": "0+"} for k in range(nz) for j in range(ny)]

#dirichlet_bc = {}
#neumann_bc = []

config = {
  "FiniteElementMethod" : {
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": [nx, ny, nz],
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "none",
    "maxIterations": 1e4,
    "bulkModulus": 1.5,
    "shearModulus": 2.0,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/out", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},
      {"format": "PythonFile", "filename": "out/out", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ]
  },
}

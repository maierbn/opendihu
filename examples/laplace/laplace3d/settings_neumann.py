import numpy as np
import sys

# 3D laplace problem

nx = 1
ny = 1
nz = 4

local = False

# Neumann boundary conditions
bc = []

# global boundary conditions
if not local:   
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)
    
      # z- plane (bottom)
      bc_value = np.sin(x*np.pi)
      bc_value = -1.0   # negative = inflow
      
      bc.append({"element": element_no, "constantValue": bc_value, "face": "2-"})
      
      # z+ plane (top)
      element_no += int((nz-1)*nx*ny)
      bc_value = np.sin(y*np.pi) + 2.0
      bc_value = 1.0
      bc.append({"element": element_no, "constantValue": bc_value, "face": "2+"})

rank_no = (int)(sys.argv[-2])

n_elements = [nx,ny,nz]

# local boundary conditions
if local:
  n_elements = [1,1,1]
  if rank_no == 0:
    n_elements = [1,1,2]

  # boundary conditions
  bc = {}
  if rank_no == 0:
    bc = {dof:1.0 for dof in range(4)}
  elif rank_no == 2:
    bc = {-1-dof:2.0 for dof in range(4)}

config = {
  "FiniteElementMethod" : {
    "nElements": n_elements,
    "nRanks": [1,1,3],
    "inputMeshIsGlobal": not local,
    "physicalExtent": [1.0, 1.0, 3.0],
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": {},
    "neumannBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "none",
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True},      
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True}
    ]
  },
}

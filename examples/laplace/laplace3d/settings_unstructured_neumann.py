import numpy as np
import sys

# 3D laplace problem

# run as (both with either local=True or local=False):
# mpirun -n 2 ./laplace ../settings_neumann.py
# ./laplace ../settings_neumann.py

local = True

nx = 2
ny = 2
nz = 5


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
n_ranks = (int)(sys.argv[-1])

n_elements = [nx,ny,nz]

# local boundary conditions
if local:
  if n_ranks > 1:
    n_elements = [2,2,2]
    if rank_no == 0:
      n_elements = [2,2,3]

  # boundary conditions
  bc = []
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)
      
      if rank_no == 0:
        bc.append({"element": element_no, "constantValue": -1.0, "face": "2-"})
        
      if rank_no == n_ranks-1:
        bc.append({"element": -(nx*ny)+element_no, "constantValue": 1.0, "face": "2+"})

# node positions
node_positions = []
for k in range(int(nz+1)):
  for j in range(int(ny+1)):
    for i in range(int(nx+1)):
      node_positions.append([i * 1.0/nx, j * 1.0/ny, k * 1.0/nz])

# elements
elements = []
for k in range(int(nz)):
  for j in range(int(ny)):
    for i in range(int(nx)):

      elements.append([k*(nx+1)*(ny+1) + j*(nx+1) + i, k*(nx+1)*(ny+1) + j*(nx+1) + i+1, k*(nx+1)*(ny+1) + (j+1)*(nx+1) + i, k*(nx+1)*(ny+1) + (j+1)*(nx+1) + i+1, \
                       (k+1)*(nx+1)*(ny+1) + j*(nx+1) + i, (k+1)*(nx+1)*(ny+1) + j*(nx+1) + i+1, (k+1)*(nx+1)*(ny+1) + (j+1)*(nx+1) + i, (k+1)*(nx+1)*(ny+1) + (j+1)*(nx+1) + i+1])

config = {
  "FiniteElementMethod" : {
    "nElements": n_elements,
    "nRanks": [1,1,n_ranks],
    "nodePositions": node_positions,
    "elements": elements,
    "inputMeshIsGlobal": not local,
    "physicalExtent": n_elements,
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": {0:0} if rank_no == 0 else {},
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

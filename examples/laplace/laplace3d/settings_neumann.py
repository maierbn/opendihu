import numpy as np
import sys

# 3D laplace problem

# run as (both with either local=True or local=False):
# mpirun -n 2 ./laplace ../settings_neumann.py
# ./laplace ../settings_neumann.py

local = True

nx = 3
ny = 3
nz = 3


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

nRanksPerCoordinateDirection = [1,1,n_ranks]

# local boundary conditions
if local:
  #if n_ranks > 1:
    #n_elements = [2,2,2]
    #if rank_no == 0:
    #  n_elements = [2,2,3]

  nRanksPerCoordinateDirection = [2,2,1]

  # boundary conditions
  bc = []
  for j in range(int(ny)):
    for i in range(int(nx)):
      x = i/nx
      y = j/ny
      element_no = int(j*nx + i)
      
      #if rank_no == 0:
      bc.append({"element": element_no, "constantValue": -1.0, "face": "2-"})
        
      #if rank_no == n_ranks-1:
      bc.append({"element": -(nx*ny)+element_no, "constantValue": 1.0, "face": "2+"})

config = {
  "FiniteElementMethod" : {
    "nElements": n_elements,
    "nRanks": nRanksPerCoordinateDirection,
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

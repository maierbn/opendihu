import numpy as np
import sys

# 3D laplace problem with quadratic Lagrange basis, run with 4 processes

# number of elements
nx = 50
ny = 50
nz = 50

# number of nodes
mx = 2*nx + 1
my = 2*ny + 1
mz = 2*nz + 1

# Dirichlet boundary conditions
bc = {}

for j in range(my):
  for i in range(mx):
    # x = 0 plane
    x = i/(2*nx)
    y = j/(2*ny)
    index = j*mx + i
  
    # z- plane (bottom)
    bc[index] = 1.0
    
    # z+ plane (top)
    index = (mz-1)*mx*my + j*mx + i
    bc[index] = 2.0
    
config = {
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "scenarioName":                   "laplace",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   None,                       # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
    # mesh parameters
    "nElements":          [nx,ny,nz],
    "nRanks":             [2,2,1],
    "inputMeshIsGlobal":  True,
    "physicalExtent":     [4.0, 4.0, 4.0],
    "outputInterval":     1.0,
    "slotName":           None,
    
    # problem parameters
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions":   [],
    "prefactor": 1,
    
    # solver parameters
    "solverType": "gmres",
    "preconditionerType": "gamg",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        

    # gamg specific options:
    "gamgType": "agg",        # one of agg, geo, or classical 
    "cycleType": "cycleV",
    "nLevels": 25,

    "maxIterations": 10000,
    "dumpFormat": "matlab",             # output solution, rhs and system matrix to a file given under dumpFilename
    "dumpFilename": "",			        		# "" means no output, insert "out" to create the output files

    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},      
      #{"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"}
    ]
  },
}

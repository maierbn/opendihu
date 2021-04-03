# 2D Laplace, quadratic elements

import numpy as np

# number of elements
nx = 10
ny = nx

# number of nodes (for quadratic elements)
mx = 2*nx+1  
my = 2*ny+1

# boundary conditions 
bc = {}
for i in range(mx):
  x = i/mx
  
  # bottom line
  bc[i] = np.sin(x*np.pi)
  
  # top line
  i2 = (my-1)*mx + i
  bc[i2] = np.sin(x*np.pi)
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "laplace",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   None,                       # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod": {
    # mesh parameters
    "nElements":                   [nx, ny],                    # number of elements in x and y direction
    "inputMeshIsGlobal":           True,                        # if the number of elements is interpreted as global number
    "physicalExtent":              [1.0, 1.0],                  # the physical (square) size of the domain
    "physicalOffset":              [0, 0],                      # the physical location of the point (0,0)
    
    # problem parameters
    "dirichletBoundaryConditions": bc,                          # the Dirichlet boundary conditions as dict
    "dirichletOutputFilename":     None,                        # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions":   [],                          # the Neumann boundary conditions as list for the elements
    "prefactor":                   1,                           # constant prefactor of the Laplace operator
    
    # solver parameters
    "solverType":                  "gmres",                     # linear solver scheme
    "preconditionerType":          "none",                      # preconditioner scheme
    "relativeTolerance":           1e-15,                       # stopping criterion, relative tolerance of the residual norm compared to inital residual
    "absoluteTolerance":           1e-10,                       # stopping criterion, absolute tolerance of the residual norm
    "maxIterations":               1e4,                         # maximum number of iterations
    "dumpFormat":                  "default",                   # format for data dump
    "dumpFilename":                None,                        # filename for dump, set to "" or None to disable
    "slotName":                    None,                        # slot name if this solver is connected to other solversr
    
    # output writers
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/laplace", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":False, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ]
  },
}

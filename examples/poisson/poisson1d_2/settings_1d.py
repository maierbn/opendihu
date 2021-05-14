import numpy as np

# Poisson 1D
n = 2  # n elements
m = n+1  # n nodes

import opendihu
# settings for quadratic discretization
if "quadratic" in opendihu.program_name:
  m = 2*n+1
  print("quadratic settings")
  
# Dirichlet boundary conditions
bc = {}
#bc[0] = 1.0
bc[0] = 0.0
bc[m-1] = 1.0

# right hand side
rhs = np.zeros(m)
for i in range(m):
  x = i / (m-1)
  rhs[i] = 6*x
  
# settings for Hermite discretization
if "hermite" in opendihu.program_name:
  m = 2*n+2
  bc = {0: 0, 1: 0, m-2: 1, m-1:3/n}  # Dirichlet bc
  #bc = {0: 0, 1:0, m-2: 1}  # Dirichlet bc
  
  # set right hand side, only nodal dofs, not derivatives
  rhs = np.zeros(m)
  
  for element_no in range(n+1):
    x = float(element_no)/n
    #rhs[2*element_no] = -10
    rhs[2*element_no] = 6*x
    rhs[2*element_no+1] = 6/n
    
print("  rhs: {}".format(rhs))

  
config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "poisson",                  # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "FiniteElementMethod" : {
    # mesh parameters
    "nElements": [n],
    "physicalExtent": [1],
    "inputMeshIsGlobal": True,
    
    # problem parameters
    "prefactor": 1,
    "rightHandSide": list(rhs),         # provide the rhs vector as list, not as numpy array
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename": None,    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    
    # solver parameters
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-15,         # 1e-10 absolute tolerance of the residual    
    "maxIterations": 10000,
    "dumpFormat": "ascii",
    "dumpFilename": "out/poisson",
    "slotName": "",
    
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out/p", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary": True, "onlyNodalValues": False, "fileNumbering": "incremental"}
    ]
  },
}

# Diffusion 2D
# number of elements
nx = 20   
ny = nx

# number of nodes
mx = 2*nx + 1
my = 2*ny + 1

# initial values
iv = {}

for j in range(int(0.2*my), int(0.3*my)):
  for i in range(int(0.5*mx), int(0.8*mx)):
    index = j*mx + i
    iv[index] = 1.0

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "diffusion",                # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "ExplicitEuler" : {
    "initialValues": iv,
    "timeStepWidth": 1e-3,
    "endTime": 0.2,
    "timeStepOutputInterval": 100,
    "checkForNanInf": False,
    "inputMeshIsGlobal": True,
    "nAdditionalFieldVariables": 0,
    "additionalSlotNames": [],
    "dirichletBoundaryConditions": {},
    "dirichletOutputFilename": None,      # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    "FiniteElementMethod" : {
      # mesh parameters
      "nElements":         [nx, ny],      # number of elements in x and y direction
      "physicalExtent":    [4.0, 4.0],    # the size of the domain in physical space
      "physicalOffset":    [0.0, 0.0],    # coordinates of lower left point
      "inputMeshIsGlobal": True,
      
      # problem parameters
      "prefactor":         0.1,           # prefactor c
      
      # solver parameters
      "solverType": "gmres",
      "preconditionerType": "none",
      "relativeTolerance": 1e-15,         # relative tolerance of the residual normal, respective to the initial residual norm, linear solver
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations": 10000,
      "dumpFormat": "default",
      "dumpFilename": "",
      "slotName": "",
    },
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "frequency": 100, "fileNumbering": "incremental"},
      {"format": "PythonFile", "outputInterval": 10, "filename": "out/diffusion2d", "binary": True, "onlyNodalValues": True, "combineFiles": True, "fileNumbering": "incremental"}
    ]
  },
}

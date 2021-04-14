# Anisotropic Diffusion 2D
n = 40   # number of elements

# initial values
iv = {}

for y in range(int(0.45*n), int(0.55*n)):
  for x in range(int(0.5*n), int(0.8*n)):
    i = y*(n+1) + x
    iv[i] = 1.0

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "diffusion",                # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "Solvers": {
    "linearSolver":
    {
      "solverType":        "gamg",
      "preconditionerType": "none",
      "relativeTolerance": 1e-15,
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
      
      # gamg specific options:
      "gamgType": "classical",        # one of agg, geo, or classical 
      "cycleType": "cycleV",
      "nLevels": 25,
      
      "maxIterations": 1000,
      "dumpFormat": "default",
      "dumpFilename": "",
    }
  },
  "ExplicitEuler" : {
    "initialValues": iv,
    "timeStepWidth": 1e-3,
    "endTime": 5.0,
    "inputMeshIsGlobal": True,
    "timeStepOutputInterval": 100,
    "checkForNanInf": False,
    "nAdditionalFieldVariables": 0,
    "additionalSlotNames": [],
    "dirichletBoundaryConditions": {},
    "dirichletOutputFilename": None,    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    "FiniteElementMethod" : {
      # mesh parameters
      "nElements":        [n,n],
      "inputMeshIsGlobal": True,
      "physicalExtent":   [4.0, 4.0],
      "physicalOffset":   [0.0, 0.0],
      
      "solverName":       "linearSolver",
      "slotName": "",
      
      # problem parameters
      "prefactor":        0.1,
      "diffusionTensor": [
        0.2, 0.2,
        0.2, 1.2
      ],
    },
    
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "frequency": 100},
      {"format": "PythonFile", "outputInterval": 100, "filename": "out/out_diffusion2d", "binary": True, "onlyNodalValues": True, "combineFiles": True, "fileNumbering": "incremental"}
    ]
  },
}

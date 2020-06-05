# Anisotropic Diffusion 2D
n = 40   # number of elements

# initial values
iv = {}

for y in range(int(0.2*n), int(0.3*n)):
  for x in range(int(0.5*n), int(0.8*n)):
    i = y*(n+1) + x
    iv[i] = 1.0

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "Solvers": {
    "linearSolver":
    {
      "solverType":        "gamg",
      "preconditionerType": "none",
      "relativeTolerance": 1e-15,
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
    }
  },
  "ExplicitEuler" : {
    "initialValues": iv,
    "timeStepWidth": 1e-3,
    "endTime": 5.0,
    
    "FiniteElementMethod" : {
      "nElements":        [n,n],
      "physicalExtent":   [4.0,4.0],
      "solverName":       "linearSolver",
      "prefactor":        0.1,
      "diffusionTensor": [
        0.2, 0.2,
        0.2, 1.2
      ],
    },
    
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "frequency": 100},
      {"format": "PythonFile", "outputInterval": 100, "filename": "out/out_diffusion2d", "binary": False}
    ]
  },
}

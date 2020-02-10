# Diffusion 2D
n = 20   # number of elements

# initial values
iv = {}

for y in range(int(0.2*n), int(0.3*n)):
  for x in range(int(0.5*n), int(0.8*n)):
    i = y*(2*n+1) + x
    iv[i] = 1.0

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "ExplicitEuler" : {
    "initialValues": iv,
    "timeStepWidth": 1e-3,
    "endTime": 1.0,
    
    "FiniteElementMethod" : {
      "nElements": [n,n],
      "physicalExtent": [4.0,4.0],
      "relativeTolerance": 1e-15,
      "prefactor": 0.1,
    },
    "OutputWriter" : [
      #{"format": "Paraview", "interval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "frequency": 100},
      {"format": "PythonFile", "outputInterval": 10, "filename": "out/diffusion2d", "binary": True}
    ]
  },
}

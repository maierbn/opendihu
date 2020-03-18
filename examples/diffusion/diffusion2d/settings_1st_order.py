# Diffusion 2D
n = 20   # number of elements

# initial values
iv = {}

for y in range(int(0.2*n), int(0.3*n)):
  for x in range(int(0.5*n), int(0.8*n)):
    i = y*(n+1) + x
    iv[i] = 1.0

print("iv: ",iv)

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "Heun" : {
    "initialValues": iv,
    "timeStepWidth": 1e-3,
    "endTime": 1.0,
    
    "FiniteElementMethod" : {
      "nElements": [n,n],
      "physicalExtent": [4.0,4.0],
      "relativeTolerance": 1e-15,    
      "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
      "dumpFilename": "out/",
      "dumpFormat": "ascii",  # ascii, default or matlab
      "prefactor": 0.1,
    },
    
    "OutputWriter" : [
      {"format": "Paraview", "interval": 1, "filename": "out/filename", "binary": "false", "fixedFormat": False, "frequency": 10},
      {"format": "PythonFile", "outputInterval": 10, "filename": "out/out_diffusion2d", "binary": True}
    ]
  },
}

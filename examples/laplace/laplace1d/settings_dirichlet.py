from datetime import datetime

# Laplace 1D Hermite

n = 5

# Dirichlet boundary conditions
bc = {}
bc[0] = 1.0
bc[1] = -1./5
bc[2*n] = 0.0
bc[2*n+1] = -1./5

config = {
  "logFormat":  "csv", # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "meta": {            # Under the "meta" key you can add any values, they will appear in the log file
    "key0": datetime.now(),
    "key1": "test",
    "key2": [0,1,2],
    "key3": {"a": "b", 1: "c"},
    "key4": None,
    "key5": 3.1415926,
  },
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 5.0,
    "dirichletBoundaryConditions": bc,
    "nodePositions": [[float(i)/n] for i in range(n+1)],
    "elements": [[i, i+1] for i in range(n)],
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "inputMeshIsGlobal": True,
    "prefactor": 1.0,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "combineFiles": False},
      {"format": "PythonFile", "outputInterval": 1, "filename": "out/p", "combineFiles": False, "onlyNodalValues": False}
    ]
  },
}

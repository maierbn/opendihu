# Diffusion 1D

# args: <n_elements> <solver> <preconditioner> <scenario_name>
#end_time = 100.0

#import numpy as np
#import pickle
import sys

#parse arguments
n = 100
n = (int)(sys.argv[0]) # number of elements
dt=0.01
end_time=1.0
solver_type = "gmres"
solver_type = sys.argv[1]
preconditioner_type = "none"
preconditioner_type = sys.argv[2]

scenario_name = sys.argv[3]
config = {
  "scenarioName": scenario_name,
  "ImplicitEuler" : {
    "preconditionerType": preconditioner_type,
    "solverType": solver_type,
    "nLevels": 25,
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 1000,
    "endTime": 1,
    "FiniteElementMethod" : {
      "preconditionerType": preconditioner_type,
      "solverType": solver_type,
      "nLevels": 25,
      #"hypreOptions": hypreOptions,
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/diffusion1d", "outputInterval": 10, "binary":False}
    ]
  },
}

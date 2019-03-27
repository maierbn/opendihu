# Diffusion 1D
n = 10000 # number of elements
dt=0.01
end_time=1
hypreOptions = "-pc_hypre_type boomeramg"

config = {
  "CrankNicolson" : {
    "preconditionerType": "none",
    "solverType": "gmres",
    "nLevels": 25,
    "hypreOptions": hypreOptions,
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 1000,
    "endTime": 1,
    "FiniteElementMethod" : {
      "preconditionerType": "none",
      "solverType": "gmres",
      "nLevels": 25,
      "hypreOptions": hypreOptions,
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

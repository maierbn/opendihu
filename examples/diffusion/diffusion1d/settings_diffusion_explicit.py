# Diffusion 1D
n = 5   # number of elements

config = {
  "ExplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 0.1,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/diffusion1d", "outputInterval": 1, "binary":False}
    ]
  },
}

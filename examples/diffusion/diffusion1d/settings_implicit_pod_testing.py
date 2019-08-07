# Diffusion 1D
n = 5   # number of elements
k = 5   

config = {
  "ModelOrderReduction": {
    "nRowsSnapshots" : n,
    "nReducedBases" : k,   
    "ImplicitEuler" : {
      "numberTimeSteps": 5,
      "endTime": 0.1,
      "initialValues": [2,2,4,5,2,2],      
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
    "ImplicitEulerReduced" : {
      "numberTimeSteps": 5,
      "endTime": 0.1,
      "initialValues": [2,2,4,5,2,2],      
      "FiniteElementMethod" : {
        "nElements": n,
        "physicalExtent": 4.0,
        "relativeTolerance": 1e-15,
        "diffusionTensor": [5.0],
      },
      "OutputWriter" : [
         #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
        {"format": "PythonFile", "filename": "out/diffusion1d_pod", "outputInterval": 1, "binary":False}
      ]
    },
  },
}

# Diffusion 1D
n = 20   # number of elements
k = 20   

config = {
  "ModelOrderReduction": {
    "nRowsSnapshots" : n,
    "nReducedBases" : k, 
    "ExplicitEuler" : {
      "numberTimeSteps": 1,
      "endTime": 1e-3,
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
    "ExplicitEulerReduced" : {
      "numberTimeSteps": 1,
      "endTime": 1e-3,
      "initialValues": [2,2,4,5,2,2],      
      "FiniteElementMethod" : {
        "nElements": n,
        "physicalExtent": 4.0,
        "relativeTolerance": 1e-15,
        "diffusionTensor": [5.0],
      },
      "OutputWriter" : [
         #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
        {"format": "PythonFile", "filename": "out/diffusion1dReduced", "outputInterval": 1, "binary":False}
      ]
    },
  },
}

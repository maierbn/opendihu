# Diffusion 1D
n = 20   # number of elements
dt=0.01
end_time=1;

config = {
  "multigrid_Vcycle": {
  "timeStepWidth": 1e-1
    "endTime": end_time,
    "durationLogKey": "duration_total",
    "timeStepOutputInterval": 10,
  "Term1": {
  "ExplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "timeStepWidth":dt/2,
    "FiniteElementMethod" : {
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
  },
  "Term2": {  
  "ExplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "timeStepWidth":dt/2,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out/diffusion1d", "outputInterval": 10, "binary":False}
    ]
  }
  }
    }
}

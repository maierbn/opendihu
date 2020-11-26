# Diffusion 1D POD
n = 5   # number of elements
k = 5   

config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "diffusion",                # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "ModelOrderReduction": {
    "nRowsSnapshots" : n,
    "nReducedBases" : k,
    "snapshots" :"./out_snapshots/snapshots.csv",
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
         {"format": "PythonFile", "filename": "out/diffusion1d_pod_full", "outputInterval": 1, "binary":False}
       ]
    },   
    "ImplicitEulerReduced" : {
      "numberTimeSteps": 5,
      "endTime": 0.1,
      "initialValues": [2,2,4,5,2,2], 
      "checkForNanInf": True,     
      "FiniteElementMethod" : {
        "nElements": n,
        "physicalExtent": 4.0,
        "relativeTolerance": 1e-15,
        "diffusionTensor": [5.0],
      },
      "OutputWriter" : [
        #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False, "fileNumbering": "incremental"},
        {"format": "PythonFile", "filename": "out/diffusion1d_pod_reduced", "outputInterval": 1, "binary":False, "fileNumbering": "incremental"}
      ]
    },
  },
}

# CellML debug
    
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def setParameters(n_instances, time_step_no, current_time, parameters, *additional_parameters):
  
  center_node = int(n_instances / 2)
  if current_time - int(current_time) < 0.1:
    print("set I_Stim for node {} to 400".format(center_node))
    parameters[center_node] = 400.
    
xdata = []
gamma_data = []
vm_data = []
      
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "ExplicitEuler" : {
    "timeStepWidth": 1e-2,
    "endTime" : 100.0,
    "initialValues": [],
    "timeStepOutputInterval": 2.5e1,
    
    "OutputWriter" : [
      #{"format": "Callback", "outputInterval": 1e4, "callback": callback},
      #{"format": "Paraview", "filename": "out", "binaryOutput": "false", "fixedFormat": False, "outputInterval": 1},
      {"format": "PythonFile", "filename": "out/cell_e", "outputInterval": 2.5e1, "binary": True}
    ],

    "CellML" : {
      "sourceFilename": "cellmlcode.cpp",
      "simdSourceFilename" : "simdcode.cpp",
      "libraryFilename": "cellml_simd_lib.so",
      "setParametersFunction": setParameters,
      "setParametersCallInterval": 1,
       #"handleResultFunction": handleResult,
      "handleResultCallInterval": 1e4,
      
      "numberStates": 57,
      "numberIntermediates": 1,   # intermediates=wanted: gamma
      "numberParameters": 2,      # parameters=known: I_Stim, l_hs
      "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
      "parametersInitialValues": [400.0, 1.0],      # parameters: I_Stim, l_hs
    },
  },
}

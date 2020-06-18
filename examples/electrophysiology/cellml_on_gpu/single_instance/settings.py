# CellML debug
    
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def setParameters(n_instances, time_step_no, current_time, parameters, *additional_parameters):
  
  center_node = int(n_instances / 2)
  if current_time - int(current_time) < 0.1 and current_time <= 3.1:
    print("set I_Stim for node {} to 400".format(center_node))
    parameters[center_node] = 400.
    
xdata = []
gamma_data = []
vm_data = []

It_value=2
tsw=0.01*pow(2,2-It_value)
et=0.3#1.0
opiv=1*pow(2,It_value-2)#10*pow..
hrciv=opiv

config = {
  "logFormat": "csv",
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "ExplicitEuler" : {
    "timeStepWidth": tsw,
    "endTime" : et,#1.0, # 0.1
cde    "initialValues": [],
    "timeStepOutputInterval": 50e1, # ignored, see line 31 instead. --aaron
    
    "OutputWriter" : [
      #{"format": "Callback", "outputInterval": 1e4, "callback": callback},
      #{"format": "Paraview", "filename": "out", "binary": "false", "fixedFormat": False, "outputInterval": 1},
      {"format": "PythonFile", "filename": "out/N_1/cell", "outputInterval": opiv, "binary": False} #0.01e1 
    ],

    "CellML" : {
      #"modelFilename": "cellmlcode.cpp",
      "gpuSourceFilename": "gpucode.cpp",
      #"simdSourceFilename" : "simdcode.cpp",
      "libraryFilename": "cellml_gpu_lib.so",
      "setParametersFunction": setParameters,
      "setParametersCallInterval": 1,
       #"handleResultFunction": handleResult,
      "handleResultCallInterval": hrciv,
      
      "statesInitialValues": [ -79.974, -80.2, 5.9, 150.9, 5.9, 12.7, 132.0, 133.0, 0.009466, 0.9952, 0.0358, 0.4981, 0.581, 0.009466, 0.9952, 0.0358, 0.4981, 0.581, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 1500.0, 0.1, 1500.0, 25, 615.0, 615.0, 811.0, 811.0, 16900.0, 16900.0, 0.4, 0.4, 7200.0, 7200.0, 799.6, 799.6, 1000.0, 1000.0, 3.0, 0.8, 1.2, 3.0, 0.3, 0.23, 0.23, 0.23, 0.23] ,

      "numberStates": 57,
      "numberAlgebraics": 1,   # algebraics=wanted: gamma
      "numberParameters": 2,      # parameters=known: I_Stim, l_hs
      "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
      "parametersInitialValues": [400.0, 1.0],      # parameters: I_Stim, l_hs
    },
  },
}

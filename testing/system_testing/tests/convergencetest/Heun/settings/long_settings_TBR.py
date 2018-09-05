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

It_value=TBR
tsw=0.001*pow(2,2-It_value)
et=100.0
opiv=1000*pow(2,It_value-2)
hrciv=opiv

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "Heun" : {
    "timeStepWidth": tsw,
    "endTime" : et,#1.0,
    "initialValues": [],
    "timeStepOutputInterval": 50eTBR, # ignored, see line 31 instead. --aaron
    
    "OutputWriter" : [
      #{"format": "Callback", "outputInterval": 1e4, "callback": callback},
      #{"format": "Paraview", "filename": "out", "binary": "false", "fixedFormat": False, "outputInterval": 1},
      {"format": "PythonFile", "filename": "out/N_TBR/cell", "outputInterval": opiv, "binary": True}
    ],

    "CellML" : {
      "sourceFilename": "cellmlcode.cpp",
      "simdSourceFilename" : "simdcode.cpp",
      "libraryFilename": "cellml_simd_lib.so",
      "setParametersFunction": setParameters,
      "setParametersCallInterval": 1,
       #"handleResultFunction": handleResult,
      "handleResultCallInterval": hrciv,
      
      "numberStates": 57,
      "numberIntermediates": 1,   # intermediates=wanted: gamma
      "numberParameters": 2,      # parameters=known: I_Stim, l_hs
      "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
      "parametersInitialValues": [400.0, 1.0],      # parameters: I_Stim, l_hs
    },
  },
}

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

It_value=4
tsw=0.001*pow(2,1-It_value)
et=0.001#1.0
opiv=1*pow(2,It_value-2)
hrciv=opiv
ElemNo=tobereplaced-1
# 7700 als kleinste Einhiet. --> 7700 *2^i --> i=9:  3942400
# ElemNo maximal ca  3.946.250 auf device 0
# ElemNo maximal ca  7.892.500 auf device 1 (doppelt so viel)

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "ExplicitEuler" : {
    #"gpu": "dummy",
    "timeStepWidth": tsw,
    "endTime" : et,#1.0, # 0.1
    "initialValues": [],
    "timeStepOutputInterval": 50e1, # ignored, see line 31 instead. --aaron
    #"logTimeStepWidthAsKey": "dt",
    "durationLogKey": "duration_timestepping",
    #"logFileName": "logs/log"+str(ElemNo),
    
    #"OutputWriter" : [
      #{"format": "Callback", "outputInterval": 1e4, "callback": callback},
      #{"format": "Paraview", "filename": "out", "binary": "false", "fixedFormat": False, "outputInterval": 1},
      #{"format": "PythonFile", "filename": "out/N_1/cell", "outputInterval": opiv, "binary": False} #0.01e1 
    #],

    "CellML" : {
      "deviceNumber": "0",
      "openaccClause": "kernels",
      "sourceFilename": "cellml_rhs.c",
      # use "gpuSourceFilename" if you want gpu offloading
      "compilerFlags": "-fPIC -ta=host,tesla:managed:cc35,cc60,time,cuda10.0 -shared -acc -I/usr/local/home/kraemer/offloading/pgi_gcc7.2.0/linux86-64/2018/mpi/openmpi-2.1.2/include",# -Minfo=accel",
      "gpuSourceFilename": "gpucodekernels.c",
      # use "libraryFilename" if you want to use existing library, needs ("useGivenLibrary": True,) as well
      #"libraryFilename": "lib/cellml_rhs_101.so", 
      #"useGivenLibrary":True,
      "setParametersFunction": setParameters,
      "setParametersCallInterval": 1,
       #"handleResultFunction": handleResult,
      "handleResultCallInterval": hrciv,
      "nElements": ElemNo,

      "statesInitialValues": [ -79.974, -80.2, 5.9, 150.9, 5.9, 12.7, 132.0, 133.0, 0.009466, 0.9952, 0.0358, 0.4981, 0.581, 0.009466, 0.9952, 0.0358, 0.4981, 0.581, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 1500.0, 0.1, 1500.0, 25, 615.0, 615.0, 811.0, 811.0, 16900.0, 16900.0, 0.4, 0.4, 7200.0, 7200.0, 799.6, 799.6, 1000.0, 1000.0, 3.0, 0.8, 1.2, 3.0, 0.3, 0.23, 0.23, 0.23, 0.23] ,
      
      "numberStates": 57,
      "numberIntermediates": 1,   # intermediates=wanted: gamma
      "numberParameters": 2,      # parameters=known: I_Stim, l_hs
      "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
      "parametersInitialValues": [400.0, 1.0],      # parameters: I_Stim, l_hs
    },
  },
}

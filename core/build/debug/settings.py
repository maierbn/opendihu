# Electrophysiology debug
n = 5

#import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58           # membrane capacitance [uF/cm^2]


print ("prefactor: ",Conductivity/(Am*Cm))

fig = plt.figure(1)
#plt.ion()

def setParameters(n_instances, time_step_no, current_time, parameters):
  #print "       > called at ",time_step_no,", time=",current_time, ", p=",parameters[0]
  
  center_node = n_instances / 2
  center_node = 0
  
  parameters[0] = 0
  if current_time - int(current_time) < 0.1 and current_time < 10:
    print ("parameters len: {}, set I_Stim for node {} to 1200".format(len(parameters),center_node))
    parameters[center_node] = 1200.
    

def handleResult(nInstances, timeStepNo, currentTime, states, intermediates):
  print ("handleResult: time step {}, t={}, nInstances: {}, n entries states: {}".format(timeStepNo, currentTime, nInstances, len(states)))
  
  xdata = []
  vm_data = []
  gamma_data = []
  for i in range(nInstances):
    xdata.append(i)
    vm_data.append(states[i])
    gamma_data.append(intermediates[i])
  
  plt.figure(1)
  plt.clf()
  plt.xlabel('position $x$')
  ax1 = plt.gca()
  ax1.plot(xdata, vm_data, 'go-', label='$V_m$')
  plt.xlabel('t')
  plt.ylabel('$V_m$')
  ax2 = ax1.twinx()
  ax2.plot(xdata, gamma_data, 'ro-', label='$\gamma$')
  plt.ylabel('$\gamma$')    
  
  # ask matplotlib for the plotted objects and their labels
  lines, labels = ax1.get_legend_handles_labels()
  lines2, labels2 = ax2.get_legend_handles_labels()
  ax2.legend(lines + lines2, labels + labels2, loc=0)
  
  filename = "out_"+str(round(currentTime,1))+".png"
  plt.savefig(filename)
  print ("   saved ""{}""".format(filename))
  #plt.draw()
    
def callback(data, shape, nEntries, dim, timeStepNo, currentTime):
  pass
    
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "Meshes": {
    "MeshFibre": {
      "nElements": 1,
      "physicalExtend": 0.1,
    },
  },
  "GodunovSplitting": {
    "disablePrinting": False,
    "disableMatrixPrinting": False,
    #"numberTimeSteps": 1,
    "timeStepWidth": 1e-1,
    "endTime": 1.0,
    "outputData1": True,
    "outputData2": False,

    "OutputWriter" : [
       #{"format": "Callback", "outputInterval": 1e4, "callback": callback},
      #{"format": "Paraview", "filename": "out", "binaryOutput": "false", "fixedFormat": False, "outputInterval": 1},
      {"format": "PythonFile", "filename": "vm", "outputInterval": 1}
    ],
    "Term1": {      # CellML
      "ExplicitEuler" : {
        "timeStepWidth": 5e-5,
        "initialValues": [],
        "timeStepOutputInterval": 1e4,
        
        "OutputWriter" : [
           #{"format": "Callback", "outputInterval": 1e4, "callback": callback},
          #{"format": "Paraview", "filename": "out", "binaryOutput": "false", "fixedFormat": False, "outputInterval": 1},
          #{"format": "PythonFile", "filename": "vm", "outputInterval": 5e4}
        ],
  
        "CellML" : {
          "sourceFilename": "cellmlcode.cpp",
          #"simdSourceFilename" : "simdcode.cpp",
          "libraryFilename": "cellml_simd_lib.so",
          #"statesInitialValues": [],
          "setParametersFunction": setParameters,
          "setParametersCallInterval": 1e3,
          "handleResultFunction": handleResult,
          "handleResultCallInterval": 2e3,
          
          "numberStates": 57,
          "numberIntermediates": 1,   # intermediates: gamma
          "numberParameters": 2,      # parameters: I_Stim, l_hs
          "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
          "parametersInitialValues": [0.0, 1.0],      # parameters: I_Stim, l_hs
          "meshName": "MeshFibre",
          "prefactor": 1.0,
        },
      },
    },
    "Term2": {     # Diffusion
      "ExplicitEuler" : {
        #"initialValues": [2,2,4,5,2,2],
        #"numberTimeSteps": 1,
        "timeStepWidth": 5e-5,
        "timeStepOutputInterval": 1e4,
        "FiniteElementMethod" : {
          #"nElements": 0,
          "physicalExtend": 1.0,
          "relativeTolerance": 1e-15,
          "meshName": "MeshFibre",
          "prefactor": Conductivity/(Am*Cm),
        },
      },
    },
  }
}

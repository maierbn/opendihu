# Electrophysiology debug
endTime = 0.1

import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58           # membrane capacitance [uF/cm^2]

n_fibers = 5

print("prefactor: ",Conductivity/(Am*Cm))

def setParameters(n_nodes, time_step_no, current_time, parameters, fibre_no):
  print("       > called setParameters at ",time_step_no,", time=",current_time, ", n_nodes=", n_nodes, ", p=",parameters[0], ", fibre ",fibre_no)
  
  return
  
  center_node = int(n_nodes / 2)
  
  parameters[0] = 0
  if current_time - int(current_time) < 0.1 and current_time < 10:
    print("parameters len: {}, set I_Stim for node {} to 1200".format(len(parameters),center_node))
    if center_node > 0:
      parameters[center_node-1] = 400.
    parameters[center_node] = 400.
    if center_node < n_nodes-1:
      parameters[center_node+1] = 400.
    
fig = plt.figure(1)
#plt.ion()

def handleResult(n_nodes, time_step_no, current_time, states, intermediates):
  #print "handleResult: time step {}, t={}, n_nodes: {}, n entries states: {}".format(time_step_no, current_time, n_nodes, len(states))
  #print "states:", states[0:n_nodes]
  
  xdata = []
  vm_data = []
  gamma_data = []
  for i in range(n_nodes):
    xdata.append(i)
    vm_data.append(states[i])
    gamma_data.append(intermediates[i])
  
  # write out a png file
  plt.figure(1)
  plt.clf()
  plt.xlabel('position $x$')
  ax1 = plt.gca()
  ax1.plot(xdata, vm_data, 'go-', label='$V_m$')
  plt.ylim(-80, 80)
  plt.xlabel('t')
  plt.ylabel('$V_m$')
  ax2 = ax1.twinx()
  ax2.plot(xdata, gamma_data, 'ro-', label='$\gamma$')
  plt.ylabel('$\gamma$')    
  plt.ylim(0, 1)
  
  # ask matplotlib for the plotted objects and their labels
  lines, labels = ax1.get_legend_handles_labels()
  lines2, labels2 = ax2.get_legend_handles_labels()
  ax2.legend(lines + lines2, labels + labels2, loc=0)
  
  filename = "out_{:06.1f}.png".format(currentTime)
  plt.savefig(filename)
  #print "   saved ""{}""".format(filename)
  #plt.draw()
    
def callback(data, shape, nEntries, dim, timeStepNo, currentTime):
  pass
    
def get_instance_config(i):
  
  instance_config = {
    "GodunovSplitting": {
      #"numberTimeSteps": 1,
      "timeStepWidth": 1e-1,
      "endTime": endTime,
      "outputData1": False,
      "outputData2": False,

      "Term1": {      # CellML
        "ExplicitEuler" : {
          "timeStepWidth": 5e-5,
          "initialValues": [],
          "timeStepOutputInterval": 1e4,
          
          "CellML" : {
            "sourceFilename": "cellmlcode.cpp",
            #"simdSourceFilename" : "simdcode.cpp",
            #"libraryFilename": "cellml_simd_lib_{}.so".format(i),
            #"statesInitialValues": [],
            "setParametersFunction": setParameters,
            "setParametersCallInterval": 1e3,
            "setParametersFunctionAdditionalParameter": i,
            #"handleResultFunction": handleResult,
            #"handleResultCallInterval": 2e3,
            #"handleResultFunctionAdditionalParameter",
            
            "numberStates": 57,
            "numberIntermediates": 1,   # intermediates: gamma
            "numberParameters": 2,      # parameters: I_Stim, l_hs
            "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
            "parametersInitialValues": [0.0, 1.0],      # parameters: I_Stim, l_hs
            "meshName": "MeshFibre"+str(i),
            "prefactor": 1.0,
          },
        },
      },
      "Term2": {     # Diffusion
        "ExplicitEuler" : {
          #"initialValues": [2,2,4,5,2,2],
          #"numberTimeSteps": 1,
          "timeStepWidth": 1e-5,
          "timeStepOutputInterval": 1e4,
          "FiniteElementMethod" : {
            #"nElements": 0,
            "physicalExtend": 1.0,
            "relativeTolerance": 1e-15,
            "meshName": "MeshFibre"+str(i),
            "prefactor": Conductivity/(Am*Cm),
          },
        },
      },
    }
  }
  return instance_config
    
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "Meshes": {
    "MeshFibre0": {
      "nElements": 2,
      "physicalExtent": 2e-3,
    },
    "MeshFibre1": {
      "nElements": 2,
      "physicalExtent": 3e-3,
    },
    "MeshFibre2": {
      "nElements": 2,
      "physicalExtent": 4e-3,
    },
  },
  "MultipleInstances": {
    "nInstances": 3,
    "instances": [get_instance_config(i) for i in range(3)],
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      {"format": "ExFile", "filename": "out/multiple_fibres", "outputInterval": 1},
      #{"format": "PythonFile", "filename": "out/"+name, "binary":False, "onlyNodalValues":True},
    ]
  }
}

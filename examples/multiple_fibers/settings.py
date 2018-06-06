# Electrophysiology debug
endTime = 10.0

import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58           # membrane capacitance [uF/cm^2]

n_fibers = 5
fibre_file = "laplace3d_structured_linear"

print("prefactor: ",Conductivity/(Am*Cm))

def setParameters(n_nodes, time_step_no, current_time, parameters, fibre_no):
  print("       > called setParameters at ",time_step_no,", time=",current_time, ", n_nodes=", n_nodes, ", p=",parameters[0], ", fibre ",fibre_no)
    
  center_node = int(n_nodes / 2)
  
  if current_time - int(current_time) < 0.1 and current_time < 10:
    print("parameters len: {}, set I_Stim for node {} to 1200".format(len(parameters),center_node))
    if center_node > 0:
      parameters[center_node-1] = 400.
    parameters[center_node] = 400.
    if center_node < n_nodes-1:
      parameters[center_node+1] = 400.
    
fig = plt.figure(1)
#plt.ion()

def debug(n_nodes, time_step_no, current_time, states, intermediates, null):
  print("time step {}, t={}, n_nodes: {}, n states: {}".format(time_step_no, current_time, n_nodes, len(states)))
  
  center_node = int(n_nodes / 2)
  
  states_of_center_node = [states[n_nodes*i + center_node] for i in range(57)]
  intermediates_of_center_node = [intermediates[n_nodes*i + center_node] for i in range(1)]
  
  print("states_of_center_node:",states_of_center_node)
  print("intermediates_of_center_node:",intermediates_of_center_node)
  
  
  node0 = int(n_nodes / 4.)
  
  states_of_node0 = [states[n_nodes*i + node0] for i in range(57)]
  intermediates_of_node0 = [intermediates[n_nodes*i + node0] for i in range(1)]
  
  print("states_of_node0:",states_of_node0)
  print("intermediates_of_node0:",intermediates_of_node0)
  wait = input("PRESS ENTER TO CONTINUE.")

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
      "timeStepWidth": 1e-3,
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
            "forceRecompileRhs": False,
            #"statesInitialValues": [],
            "setParametersFunction": setParameters,
            "setParametersCallInterval": 1e3,
            #"setParametersFunctionAdditionalParameter": i,
            
            #"handleResultFunction": debug,
            "handleResultCallInterval": 1,
            "handleResultFunctionAdditionalParameter": i,
            
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
          "timeStepWidth": 2e-7,
          "timeStepOutputInterval": 1e4,
          "FiniteElementMethod" : {
            "relativeTolerance": 1e-15,
            "meshName": "MeshFibre"+str(i),
            "prefactor": Conductivity/(Am*Cm),
          },
        },
      },
    }
  }
  return instance_config
    
    
meshes = {}

with open(fibre_file, "rb") as f:
  streamlines = pickle.load(f)
    
nInstances = len(streamlines)
print("nInstances: {}".format(nInstances))


nInstances = 1
    
for i,streamline in enumerate(streamlines):
  meshes["MeshFibre{}".format(i)] = {
    "nElements": len(streamline)-1,
    "nodePositions": streamline
  }
    
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "Meshes": meshes,
  "MultipleInstances": {
    "nInstances": nInstances,
    "instances": [get_instance_config(i) for i in range(nInstances)],
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      {"format": "ExFile", "filename": "out/multiple_fibres", "outputInterval": 1},
      #{"format": "PythonFile", "filename": "out/"+name, "binary":False, "onlyNodalValues":True},
    ]
  }
}

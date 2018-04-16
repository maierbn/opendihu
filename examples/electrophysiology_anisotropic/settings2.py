# Electrophysiology release (used for runtime measurements)

import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58           # membrane capacitance [uF/cm^2]

# adjustable parameters
nElements = 1000
endTime = 10.0
cellml_type = "non-simd"  # non-simd, native-simd, simd-opt

print "sys.argv: ",sys.argv

if len(sys.argv) > 2:
  nElements = sys.argv[2]
  
if len(sys.argv) > 3:
  cellml_type = sys.argv[3]
  
print "nElements: ",nElements,", cellml_type=",cellml_type

def setParameters(n_instances, time_step_no, current_time, parameters):
  #print "       > called at ",time_step_no,", time=",current_time, ", p=",parameters[0]
  
  center_node = n_instances / 2
  
  # set Istim to 0
  if center_node > 0:
    parameters[center_node-1] = 0.0
  parameters[center_node] = 0.0
  if center_node < n_instances-1:
    parameters[center_node+1] = 0.0
      
      
  # when condition is true, set Istim to 3x400
  if current_time - int(current_time) < 0.1 and current_time < 10:
    print "parameters len: {}, set I_Stim for node {} to 1200".format(len(parameters),center_node)
    if center_node > 0:
      parameters[center_node-1] = 400.
    parameters[center_node] = 400.
    if center_node < n_instances-1:
      parameters[center_node+1] = 400.
    
fig = plt.figure(1)
#plt.ion()

def handleResult(nInstances, timeStepNo, currentTime, states, intermediates):
  #print "handleResult: time step {}, t={}, nInstances: {}, n entries states: {}".format(timeStepNo, currentTime, nInstances, len(states))
  
  xdata = []
  vm_data = []
  gamma_data = []
  for i in range(nInstances):
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
  print "   saved ""{}""".format(filename)
  #plt.draw()
    
def callback(data, shape, nEntries, dim, timeStepNo, currentTime):
  pass
  
cellml = {
  "setParametersFunction": setParameters,
  "setParametersCallInterval": 2e2,
  "handleResultFunction": handleResult,
  "handleResultCallInterval": 2e3,
  
  "numberStates": 57,
  "numberIntermediates": 1,   # intermediates: gamma
  "numberParameters": 2,      # parameters: I_Stim, l_hs
  "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
  "parametersInitialValues": [0.0, 1.0],      # parameters: I_Stim, l_hs
  "meshName": "MeshFibre",
  "prefactor": 1.0,
}

# 1) if simdSourceFilename is given, use that source to compile the library
# 2) if not 1) but sourceFilename is given, create simdSourceFilename from that and compile library
# 3) if not 2) but libraryFilename is given, load that library, if it contains simdRhs, use that, if it contains non-simd rhs use that

if cellml_type == "non-simd":
  cellml["libraryFilename"] = "cellml_non_simd_lib.so"
elif cellml_type == "native-simd":
  cellml["sourceFilename"] = "cellmlcode.cpp"
elif cellml_type == "simd-opt":
  cellml["simdSourceFilename"] = "simd-opt.cpp"
  
    
config = {
  "Meshes": {
    "MeshFibre": {
      "nElements": nElements,
      "physicalExtend": 5.0,
    },
  },
  "GodunovSplitting": {
    "timeStepWidth": 1e-1,
    "endTime": endTime,
    "outputData1": False,
    "outputData2": False,

    "Term1": {      # CellML
      "ExplicitEuler" : {
        "timeStepWidth": 1e-5,
        "initialValues": [],
        "timeStepOutputInterval": 1e4,
  
        "CellML" : cellml,
      },
    },
    "Term2": {     # Diffusion
      "ExplicitEuler" : {
        "timeStepWidth": 5e-4,
        "timeStepOutputInterval": 1e4,
        "FiniteElementMethod" : {
          "physicalExtend": 1.0,
          "relativeTolerance": 1e-5,
          "meshName": "MeshFibre",
          "prefactor": Conductivity/(Am*Cm),
        },
      },
    },
  }
}

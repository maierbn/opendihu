# CellML release
    
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def setParameters(n_instances, time_step_no, current_time, parameters):
  #print "       > called at ",time_step_no,", time=",current_time, ", p=",parameters[0]
  
  center_node = int(n_instances / 2)
  
  parameters[0] = 0
  if current_time - int(current_time/10)*10 < 0.1 and current_time > 1000:
    #print "set I_Stim for node {} to 1200".format(center_node)
    parameters[center_node] = 1200.
    
fig = plt.figure(1)
#plt.ion()

xdata = []
gamma_data = []
vm_data = []
ca_data = []

def handleResult(nInstances, timeStepNo, currentTime, states, intermediates):
#def callback(data, shape, nEntries, dim, timeStepNo, currentTime):
  #print "time step {}, t={}, nEntries: {}, dim: {}, data: {}".format(timeStepNo, currentTime, nEntries, dim, data)
  #print "time step {}, t={}".format(timeStepNo, currentTime)
  #print "intermediates: {}".format(intermediates)
    
  if currentTime > 100:
    
    xdata.append(currentTime)
    gamma_data.append(intermediates[0])
    vm_data.append(states[0])
    ca_data.append(states[30])
    
    #print "plot xdata: {} ydata: {}".format(str(xdata), str(ydata))
    
    if False:
      plt.figure(1)
      plt.clf()
      ax1 = plt.gca()
      ax1.plot(xdata, vm_data, 'go-', label='$V_m$')
      plt.xlabel('t')
      plt.ylabel('$V_m$')
      ax2 = ax1.twinx()
      ax2.plot(xdata, gamma_data, 'ro-', label='$\gamma$')
      ax2.plot(xdata, ca_data, 'ro-', label='Ca')
      plt.ylabel('$\gamma$')    
      
      # ask matplotlib for the plotted objects and their labels
      lines, labels = ax1.get_legend_handles_labels()
      lines2, labels2 = ax2.get_legend_handles_labels()
      ax2.legend(lines + lines2, labels + labels2, loc=0)
      plt.savefig("out0.png")
    
    if int(currentTime) - int(int(currentTime)/100)*100 < 0.1:
      plt.figure(2)
      plt.clf()
      plt.subplot(311)
      plt.plot(xdata, vm_data, 'r-', label='$V_m$')
      plt.ylabel('$V_m [mV]$')
      plt.legend()
      
      plt.subplot(312)
      plt.plot(xdata, gamma_data, 'g-', label='Ca')
      plt.ylabel('Ca [-]')
      plt.legend()
      
      plt.subplot(313)
      plt.plot(xdata, ca_data, 'b-', label='$\gamma$')
      plt.ylabel('$\gamma$ [-]')
      plt.xlabel('t [ms]')
      plt.legend()
      
      plt.savefig("out.png")
    #plt.draw()
      

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "ExplicitEuler" : {
    "timeStepWidth": 1e-5,
    "endTime" : 20000.0,
    "initialValues": [],
    "timeStepOutputInterval": 1e5,
    
    "OutputWriter" : [
       #{"format": "Callback", "outputInterval": 1e4, "callback": callback},
      #{"format": "Paraview", "filename": "out", "binaryOutput": "false", "fixedFormat": False, "outputInterval": 1},
      #{"format": "Python", "filename": "p", "outputInterval": 1}
    ],

    "CellML" : {
      "sourceFilename": "cellmlcode.cpp",
      "simdSourceFilename" : "simdcode.cpp",
      "libraryFilename": "cellml_simd_lib.so",
      "setParametersFunction": setParameters,
      "setParametersCallInterval": 1e3,
      "handleResultFunction": handleResult,
      "handleResultCallInterval": 1e3,
      
      "numberStates": 57,
      "numberIntermediates": 1,   # intermediates=wanted: gamma
      "numberParameters": 2,      # parameters=known: I_Stim, l_hs
      "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma, state 28 = Ca
      "parametersInitialValues": [0.0, 1.0],      # parameters: I_Stim, l_hs
    },
  },
}

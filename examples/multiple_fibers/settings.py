# Electrophysiology debug
endTime = 100.0

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

fibre_file = "../input/laplace3d_structured_linear"
fibre_distribution_file = "../input/MU_fibre_distribution_3780.txt"
firing_times_file = "../input/MU_firing_times_real.txt"

print("prefactor: ",Conductivity/(Am*Cm))
print("numpy path: ",np.__path__)

def setParameters(n_nodes, time_step_no, current_time, parameters, fibre_no):
  
  # determine motor unit
  mu_no = fibre_distribution[fibre_no % len(fibre_distribution)]-1
  
  # determine if fibre fires now
  frequency = 10.0 # Hz
  index = int(current_time * frequency)
  n_firing_times = np.size(firing_times,0)
  fibre_gets_stimulated = firing_times[index % n_firing_times, mu_no] == 1
  
  print("       > called setParameters at timestep {}, t={}, n_nodes, IStim={}, fibre no {}, MU {}, stimulated: {}".\
    format(time_step_no, current_time, n_nodes, parameters[0], fibre_no, mu_no, fibre_gets_stimulated))
    
  if not fibre_gets_stimulated:
    pass
    #return
    
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width = 1.  # cm
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node = int(n_nodes / 2) + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate = [innervation_node]
  if innervation_node > 0:
    nodes_to_stimulate.insert(0, innervation_node-1)
  if innervation_node < n_nodes-1:
    nodes_to_stimulate.append(innervation_node)
  
  # stimulation value
  stimulation_current = 400.
  
  for node_no in nodes_to_stimulate:
    parameters[node_no] = stimulation_current
  
  print("set stimulation for nodes {}".format(nodes_to_stimulate))
  
  wait = input("Press any key to continue...")
    
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
      "outputData2": True,

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
            "setParametersCallInterval": 1,          # setParameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
            "setParametersFunctionAdditionalParameter": i,
            
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
            "relativeTolerance": 1e-10,
            "meshName": "MeshFibre"+str(i),
            "prefactor": Conductivity/(Am*Cm),
          },
          "OutputWriter" : [
            {"format": "Paraview", "outputInterval": 1e5, "filename": "out/fibre_"+str(i), "binaryOutput": True, "fixedFormat": False},
            {"format": "ExFile", "filename": "out/fibre_"+str(i), "outputInterval": 1e5},
            {"format": "PythonFile", "filename": "out/fibre_"+str(i), "outputInterval": 1e5, "binary":True, "onlyNodalValues":True},
          ]
        },
      },
    }
  }
  return instance_config
    
    
# create fibre meshes
meshes = {}

with open(fibre_file, "rb") as f:
  streamlines = pickle.load(f)
    
nInstances = len(streamlines)
print("nInstances: {}".format(nInstances))

#nInstances = 1
    
for i,streamline in enumerate(streamlines):
  meshes["MeshFibre{}".format(i)] = {
    "nElements": len(streamline)-1,
    "nodePositions": streamline
  }
    
# load MU distribution and firing times
fibre_distribution = np.genfromtxt(fibre_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

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
      {"format": "PythonFile", "filename": "out/multiple_fibres", "binary":True, "onlyNodalValues":True},
    ]
  }
}
    

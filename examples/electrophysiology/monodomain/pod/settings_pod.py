# Electrophysiology
# Monodomain with Hodgkin-Huxley model as rhs
#
# This file was created by Nehzat.
# parameters: [<n_elements> [<end_time>]]

end_time = 200.0   # [ms] end time of simulation
n_elements = 500

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 0.  # cm
cellml_file = "../input/hodgkin_huxley_1952.c"

# timing parameters
stimulation_frequency = 10.0      # [1/ms] frequency if which stimulation current can be switched on and off
dt_0D = 5e-5                      # timestep width of ODEs, cellml integration
dt_1D = 1e-5                      # timestep width of diffusion
dt_3D = 1e-1                      # overall timestep width of splitting scheme
output_timestep = 1e-1            # timestep for output files

# import needed packages
import sys
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# read number of elements from command line
if len(sys.argv) > 0:
  try:
    n_elements = int(sys.argv[0])
  except:
    print("could not parse n_elements {}".format(n_elements))
if len(sys.argv) > 1:
  try:
    end_time = float(sys.argv[1])
  except:
    print("could not parse end_time {}".format(end_time))

print("n elements: {}, end time: {}".format(n_elements,end_time))
print("prefactor: ",Conductivity/(Am*Cm))

# determine if fibre gets stimulation at given time
def fibre_gets_stimulated(current_time):
  a = current_time * stimulation_frequency
  
  if a - int(a) < 0.1 and current_time < 5:
    return True
  else:
    return False
  
# callback function that can set parameters, i.e. stimulation current
def set_parameters(n_nodes, time_step_no, current_time, parameters, null):
  
  # determine if fibre gets stimulated at the current time
  is_fibre_gets_stimulated = fibre_gets_stimulated(current_time)
  
  # determine nodes to stimulate (center node, left and right neighbour)
  innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
  innervation_node = int(n_nodes / 2) + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
  nodes_to_stimulate = [innervation_node]
  if innervation_node > 0:
    nodes_to_stimulate.insert(0, innervation_node-1)
  if innervation_node < n_nodes-1:
    nodes_to_stimulate.append(innervation_node+1)
  
  # stimulation value
  if is_fibre_gets_stimulated:
    stimulation_current = 400.
  else:
    stimulation_current = 0.
  
  for node_no in nodes_to_stimulate:
    parameters[node_no] = stimulation_current

  print("       set_parameters at timestep {}, t={}, n_nodes={}, stimulated: {}".\
    format(time_step_no, current_time, n_nodes, is_fibre_gets_stimulated))
 
  print("       set stimulation for nodes {}".format(nodes_to_stimulate))
  
  #wait = input("Press any key to continue...")
    
fig = plt.figure(1)
#plt.ion()

# callback function that is called after integration of rhs, generates plots
def handleResult(n_instances, time_step_no, current_time, states, algebraics, name_information, null):
  
  # collect data for every instance
  xdata = []
  vm_data = []
  gamma_data = []
  for i in range(nInstances):
    xdata.append(i)
    vm_data.append(states[i])
    gamma_data.append(algebraics[i])
  
  # generate plot of Vm and gamma
  # prepare figure
  plt.figure(1)
  plt.clf()
  plt.xlabel('position $x$')
  ax1 = plt.gca()
  ax1.plot(xdata, vm_data, 'go-', label='$V_m$')
  plt.ylim(-80, 80)
  plt.xlabel('t')
  plt.ylabel('$V_m$')
  ax2 = ax1.twinx()
  
  # plot data
  ax2.plot(xdata, gamma_data, 'ro-', label='$\gamma$')
  plt.ylabel('$\gamma$')    
  plt.ylim(0, 1)
  
  # ask matplotlib for the plotted objects and their labels
  lines, labels = ax1.get_legend_handles_labels()
  lines2, labels2 = ax2.get_legend_handles_labels()
  ax2.legend(lines + lines2, labels + labels2, loc=0)
  
  # save to png file
  filename = "out_{:06.1f}.png".format(currentTime)
  plt.savefig(filename)
  #print "   saved ""{}""".format(filename)
  #plt.draw()
    
# callback function from output writer
def callback(data, shape, nEntries, dim, timeStepNo, currentTime, null):
  pass
    
config = {
  "logFormat": "csv",
  "Meshes": {
    "MeshFibre": {
      "nElements": n_elements,
      "physicalExtent": n_elements/10.,
    },
  },
  "GodunovSplitting": {
    #"numberTimeSteps": 1,
    "timeStepWidth": dt_3D,  # 1e-1
    "endTime": end_time,
    "connectedSlotsTerm1To2": [0],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
    "connectedSlotsTerm2To1": [0],   # transfer the same back
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1e5, "filename": "out/fibre_splitting", "binaryOutput": True, "fixedFormat": False},
      {"format": "ExFile", "filename": "out/fibre_splitting", "outputInterval": 1e5, "sphereSize": "2*2*2"},
    ],
    "Term1": {      # CellML
      "ModelOrderReduction": {
        "nFullBasis" : n_elements,
        "nReducedBasis" : n_elements/10,
        "ExplicitEuler" : {
          "timeStepWidth": dt_0D,  # 5e-5
          "initialValues": [],
          "timeStepOutputInterval": 1e4,
  
          "CellML" : {
            "modelFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
            #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from modelFilename and is ready for multiple instances
            #"libraryFilename": "cellml_simd_lib.so",   # compiled library
            #"statesInitialValues": [],
            "setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
            "setParametersCallInterval": 1./stimulation_frequency/dt_0D,     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
            #"handleResultFunction": handleResult,
            #"handleResultCallInterval": 2e3,
          
            "statesForTransfer": 0,     # state 0 = Vm
            "parametersUsedAsAlgebraic": [],       # list of algebraic value indices, that will be set by parameters. Explicitely defined parameters that will be copied to algebraics, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersUsedAsConstant": [2],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
            "parametersInitialValues": [0.0],      # initial values for the parameters: I_Stim
            "meshName": "MeshFibre",
            "prefactor": 1.0,
          },
        
          "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": 1e4, "filename": "out/states", "binary": True},
          ],
        },
      },
    },
    "Term2": {     # Diffusion
      "ModelOrderReduction": {
        "nFullBasis" : n_elements,
        "nReducedBasis" : n_elements/10,
        "ExplicitEuler" : {
          #"initialValues": [2,2,4,5,2,2],
          #"numberTimeSteps": 1,
          "timeStepWidth": dt_1D,
          "timeStepOutputInterval": 1e4,
          "FiniteElementMethod" : {
            "physicalExtend": 1.0,
            "relativeTolerance": 1e-15,
            "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual    
            "meshName": "MeshFibre",
            "prefactor": Conductivity/(Am*Cm),
            },
            "OutputWriter" : [
              {"format": "PythonFile", "outputInterval": 1e5, "filename": "out/fibre", "binary": True},
              {"format": "Paraview", "outputInterval": 1e5, "filename": "out/fibre", "binaryOutput": True, "fixedFormat": False},
              {"format": "ExFile", "filename": "out/fibre", "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
            ],
        },
      },
    }
}

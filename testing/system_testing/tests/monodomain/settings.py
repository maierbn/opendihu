# Electrophysiology
# Monodomain with either Shorten or Hodgkin-Huxley model as rhs
#
# arguments: <name> <cellml_file> <end_time>

end_time = 1000.0   # [ms] end time of simulation
n_elements = 500

# global parameters
PMax = 7.3              # maximum stress [N/cm^2]
Conductivity = 3.828    # sigma, conductivity [mS/cm]
Am = 500.0              # surface area to volume ratio [cm^-1]
Cm = 0.58               # membrane capacitance [uF/cm^2]
innervation_zone_width = 1.  # cm
innervation_zone_width = 0.  # cm


cellml_file = "../input/shorten_ocallaghan_davidson_soboleva_2007.c"
cellml_file = "../input/shorten_opencmiss.cpp"

# timing parameters
stimulation_frequency = 10.0      # [1/ms] frequency if which stimulation current can be switched on and off
dt_0D = 5e-5                      # timestep width of ODEs, cellml integration
dt_1D = 1e-5                      # timestep width of diffusion
dt_3D = 1e-1                      # overall timestep width of splitting scheme
output_timestep = 1e-1            # timestep for output files

print("prefactor: ",Conductivity/(Am*Cm))

# import needed packages
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 2:
  if "check_results.py" not in sys.argv[0]:
    name = sys.argv[0]
    cellml_file = sys.argv[1]
    end_time = float(sys.argv[2])
    
    print("name: \"{}\", cellml_file: {}".format(name, cellml_file))


if "shorten" in cellml_file:
  parametersUsedAsIntermediate = [32]
  parametersUsedAsConstant = [65]
  parametersInitialValues = [0.0, 1.0]
  
elif "hodgkin_huxley" in cellml_file:
  parametersUsedAsIntermediate = []
  parametersUsedAsConstant = [2]
  parametersInitialValues = [0.0]
  

# determine if fibre gets stimulation at given time
def fibre_gets_stimulated(current_time):
  a = current_time * stimulation_frequency
  
  if a - int(a) < 0.1 and a < 5:
    return True
  else:
    return False
  
# callback function that can set parameters, i.e. stimulation current
def set_parameters(n_nodes, time_step_no, current_time, parameters, dof_nos_global, null):
  
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
def handleResult(n_instances, time_step_no, current_time, states, intermediates, null):
  
  # collect data for every instance
  xdata = []
  vm_data = []
  gamma_data = []
  for i in range(nInstances):
    xdata.append(i)
    vm_data.append(states[i])
    gamma_data.append(intermediates[i])
  
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
    "outputData1": False,
    "outputData2": False,

    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1e5, "filename": "out/fibre_splitting", "binary": True, "fixedFormat": False},
      {"format": "ExFile", "filename": "out/fibre_splitting", "outputInterval": 1e5, "sphereSize": "2*2*2"},
    ],
    "Term1": {      # CellML
      "ExplicitEuler" : {
        "timeStepWidth": dt_0D,  # 5e-5
        "initialValues": [],
        "timeStepOutputInterval": 1e4,
  
        "CellML" : {
          "sourceFilename": cellml_file,             # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
          #"simdSourceFilename" : "simdcode.cpp",     # transformed C++ source file that gets generated from sourceFilename and is ready for multiple instances
          #"libraryFilename": "cellml_simd_lib.so",   # compiled library
          "forceRecompileRhs": True,
          #"statesInitialValues": [],
          "setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
          "setParametersCallInterval": 1./stimulation_frequency/dt_0D,     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
          #"handleResultFunction": handleResult,
          #"handleResultCallInterval": 2e3,
          
          "outputStateIndex": 0,     # state 0 = Vm, rate 28 = gamma
          "parametersUsedAsIntermediate": parametersUsedAsIntermediate,  #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersUsedAsConstant": parametersUsedAsConstant,          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
          "parametersInitialValues": parametersInitialValues,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
          "meshName": "MeshFibre",
          "prefactor": 1.0,
        },
        
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": 1e4, "filename": "out/"+name, "binary": True},
        ],
      },
    },
    "Term2": {     # Diffusion
      "ExplicitEuler" : {
        #"initialValues": [2,2,4,5,2,2],
        #"numberTimeSteps": 1,
        "timeStepWidth": dt_1D,
        "timeStepOutputInterval": 1e4,
        "FiniteElementMethod" : {
          #"nElements": 0,
          "physicalExtend": 1.0,
          "relativeTolerance": 1e-15,
          "meshName": "MeshFibre",
          "prefactor": Conductivity/(Am*Cm),
        },
        "OutputWriter" : [
          {"format": "PythonFile", "outputInterval": 1e5, "filename": "out/vm_"+name, "binary": True},
          {"format": "Paraview", "outputInterval": 1e5, "filename": "out/vm_"+name, "binary": True, "fixedFormat": False},
          {"format": "ExFile", "filename": "out/vm_"+name, "outputInterval": 1e5, "sphereSize": "0.02*0.02*0.02"},
        ],
      },
    },
  }
}

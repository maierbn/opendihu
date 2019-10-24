#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "equation/diffusion.h"
#include "../utility.h"

TEST(CellMLTest, HodgkinHuxley)
{
  std::string pythonConfig = R"(

# timing parameters
stimulation_frequency = 10.0      # [1/ms] frequency if which stimulation current can be switched on and off
dt_0D = 5e-5                      # timestep width of ODEs, cellml integration

# CellML Hodgkin-Huxley from cpp file
config = {
  "ExplicitEuler" : {
    "timeStepWidth": 1e-5,
    "endTime" : 1.0,
    "initialValues": [],
    "timeStepOutputInterval": 1e5,

    "OutputWriter" : [
      {"format": "PythonFile", "filename": "out", "binary": False, "outputInterval": 1e4}
    ],

    "CellML" : {
      "sourceFilename": "../input/hodgkin_huxley_1952.c",
      "setParametersCallInterval": 1e3,
      "useGivenLibrary": False,
      #"statesInitialValues": [-75,  .05, 0.6, 0.325],
      "statesInitialValues": [-20, 0.05, 0.6, 0.325],
      "parametersInitialValues": [400.0],      # initial values for the parameters: I_Stim
      #"setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
      #"setParametersCallInterval": 1./stimulation_frequency/dt_0D,     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05

      "parametersUsedAsIntermediate": [],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
      "parametersUsedAsConstant": [2],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
    },
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  TimeSteppingScheme::ExplicitEuler<
    CellmlAdapter<4>
  > problem(settings);

  problem.run();

  std::string referenceOutput = "{\"meshType\": \"StructuredRegularFixed\", \"dimension\": 1, \"nElementsGlobal\": [0], \"nElementsLocal\": [0], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0]}, {\"name\": \"y\", \"values\": [0.0]}, {\"name\": \"z\", \"values\": [0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"V\", \"values\": [36.18142823585638]}, {\"name\": \"m\", \"values\": [0.9987345768519429]}, {\"name\": \"h\", \"values\": [0.2446134695357078]}, {\"name\": \"n\", \"values\": [0.5789949501440312]}]}], \"timeStepNo\": 90001, \"currentTime\": 0.90001}";
  assertFileMatchesContent("out_0000009.py", referenceOutput);
}

TEST(CellMLTest, ShortenOpenCMISS)
{
  std::string pythonConfig = R"(

# timing parameters
stimulation_frequency = 10.0      # [1/ms] frequency if which stimulation current can be switched on and off
dt_0D = 5e-5                      # timestep width of ODEs, cellml integration

# CellML Shorten from OpenCMISS generated cpp file
config = {
  "ExplicitEuler" : {
    "timeStepWidth": 1e-5,
    "endTime" : 1.0,
    "initialValues": [],
    "timeStepOutputInterval": 1e5,

    "OutputWriter" : [
      {"format": "PythonFile", "filename": "out", "binary": False, "outputInterval": 1e4}
    ],

    "CellML" : {
      "sourceFilename": "../input/shorten_opencmiss.cpp",
      "setParametersCallInterval": 1e3,
      "useGivenLibrary": False,
      "parametersUsedAsIntermediate": [32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
      "parametersUsedAsConstant": [65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
      "parametersInitialValues": [1000.0, 1.0],      # initial values for the parameters: I_Stim, l_hs

      #"setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
      #"setParametersCallInterval": 1./stimulation_frequency/dt_0D,     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
},
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  TimeSteppingScheme::ExplicitEuler<
    CellmlAdapter<57>
  > problem(settings);

  problem.run();

  std::string referenceOutput = "{\"meshType\": \"StructuredRegularFixed\", \"dimension\": 1, \"nElementsGlobal\": [0], \"nElementsLocal\": [0], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0]}, {\"name\": \"y\", \"values\": [0.0]}, {\"name\": \"z\", \"values\": [0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"vS\", \"values\": [134.25831982418282]}, {\"name\": \"vT\", \"values\": [129.46627066359946]}, {\"name\": \"K_t\", \"values\": [5.904060365578467]}, {\"name\": \"K_i\", \"values\": [150.89863228412054]}, {\"name\": \"K_e\", \"values\": [5.906803266826562]}, {\"name\": \"Na_i\", \"values\": [12.699713043228583]}, {\"name\": \"Na_t\", \"values\": [131.99288848618707]}, {\"name\": \"Na_e\", \"values\": [133.0014881414964]}, {\"name\": \"n\", \"values\": [0.8325126231102002]}, {\"name\": \"h_K\", \"values\": [0.45117680867638915]}, {\"name\": \"m\", \"values\": [0.9999988911096568]}, {\"name\": \"h\", \"values\": [0.010757587063453006]}, {\"name\": \"S\", \"values\": [0.5797670629801744]}, {\"name\": \"n_t\", \"values\": [0.6835959680178263]}, {\"name\": \"h_K_t\", \"values\": [1.54049441227643e-10]}, {\"name\": \"m_t\", \"values\": [0.9999983751266015]}, {\"name\": \"h_t\", \"values\": [0.017298326243061504]}, {\"name\": \"S_t\", \"values\": [0.5803050265994951]}, {\"name\": \"O_0\", \"values\": [4.57410818823492e-10]}, {\"name\": \"O_1\", \"values\": [2.524237133190971e-07]}, {\"name\": \"O_2\", \"values\": [3.935788241470067e-05]}, {\"name\": \"O_3\", \"values\": [0.002325163091925736]}, {\"name\": \"O_4\", \"values\": [0.0913566128987036]}, {\"name\": \"C_0\", \"values\": [0.0002519254215580167]}, {\"name\": \"C_1\", \"values\": [0.006987590727294847]}, {\"name\": \"C_2\", \"values\": [0.07235721571073575]}, {\"name\": \"C_3\", \"values\": [0.3255351150464797]}, {\"name\": \"C_4\", \"values\": [0.501146766339743]}, {\"name\": \"dummy\", \"values\": [-0.008905857405137775]}, {\"name\": \"Ca_1\", \"values\": [51.96380790766034]}, {\"name\": \"Ca_SR1\", \"values\": [60.741871506753924]}, {\"name\": \"Ca_2\", \"values\": [9.572137104634447]}, {\"name\": \"Ca_SR2\", \"values\": [123.05630923870143]}, {\"name\": \"Ca_T_2\", \"values\": [38.010886237108934]}, {\"name\": \"Ca_P1\", \"values\": [615.0]}, {\"name\": \"Ca_P2\", \"values\": [615.0]}, {\"name\": \"Mg_P1\", \"values\": [811.0]}, {\"name\": \"Mg_P2\", \"values\": [811.0]}, {\"name\": \"Ca_Cs1\", \"values\": [16879.144075133656]}, {\"name\": \"Ca_Cs2\", \"values\": [16883.28569569442]}, {\"name\": \"Ca_ATP1\", \"values\": [95.1029020400325]}, {\"name\": \"Ca_ATP2\", \"values\": [36.866205355349955]}, {\"name\": \"Mg_ATP1\", \"values\": [7227.6969156916]}, {\"name\": \"Mg_ATP2\", \"values\": [7229.597181245323]}, {\"name\": \"ATP1\", \"values\": [677.2001822683987]}, {\"name\": \"ATP2\", \"values\": [733.5366133993373]}, {\"name\": \"Mg1\", \"values\": [970.8347934879887]}, {\"name\": \"Mg2\", \"values\": [970.4176499751273]}, {\"name\": \"Ca_CaT2\", \"values\": [9.545530677855714]}, {\"name\": \"D_0\", \"values\": [0.6496734048953928]}, {\"name\": \"D_1\", \"values\": [1.126353189326496]}, {\"name\": \"D_2\", \"values\": [3.3311608670301185]}, {\"name\": \"A_1\", \"values\": [0.3109211536150155]}, {\"name\": \"A_2\", \"values\": [0.23081380116676098]}, {\"name\": \"P\", \"values\": [0.23000829237218517]}, {\"name\": \"P_SR\", \"values\": [0.23011860091881098]}, {\"name\": \"P_C_SR\", \"values\": [0.22988139450365325]}]}], \"timeStepNo\": 90001, \"currentTime\": 0.90001}";
  assertFileMatchesContent("out_0000009.py", referenceOutput);
}

TEST(CellMLTest, ShortenOpenCOR)
{
  std::string pythonConfig = R"(

# timing parameters
stimulation_frequency = 10.0      # [1/ms] frequency if which stimulation current can be switched on and off
dt_0D = 5e-5                      # timestep width of ODEs, cellml integration

# CellML Shorten from OpenCOR generated cpp file
config = {
  "ExplicitEuler" : {
    "timeStepWidth": 1e-5,
    "endTime" : 1.0,
    "initialValues": [],
    "timeStepOutputInterval": 1e5,

    "OutputWriter" : [
      {"format": "PythonFile", "filename": "out", "binary": False, "outputInterval": 1e4}
    ],

    "CellML" : {
      "sourceFilename": "../input/shorten_ocallaghan_davidson_soboleva_2007.c",
      "setParametersCallInterval": 1e3,
      "useGivenLibrary": False,
      "parametersUsedAsIntermediate": [32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
      "parametersUsedAsConstant": [65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
      "parametersInitialValues": [1000.0, 1.0],      # initial values for the parameters: I_Stim, l_hs

      #"setParametersFunction": set_parameters,    # callback function that sets parameters like stimulation current
      #"setParametersCallInterval": 1./stimulation_frequency/dt_0D,     # set_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
},
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  TimeSteppingScheme::ExplicitEuler<
    CellmlAdapter<57, 71>  // 57 states, 71 intermediates
  > problem(settings);

  problem.run();

  std::string referenceOutput = "{\"meshType\": \"StructuredRegularFixed\", \"dimension\": 1, \"nElementsGlobal\": [0], \"nElementsLocal\": [0], \"beginNodeGlobalNatural\": [0], \"hasFullNumberOfNodes\": [true], \"basisFunction\": \"Lagrange\", \"basisOrder\": 1, \"onlyNodalValues\": true, \"nRanks\": 1, \"ownRankNo\": 0, \"data\": [{\"name\": \"geometry\", \"components\": [{\"name\": \"x\", \"values\": [0.0]}, {\"name\": \"y\", \"values\": [0.0]}, {\"name\": \"z\", \"values\": [0.0]}]}, {\"name\": \"solution\", \"components\": [{\"name\": \"vS\", \"values\": [23.814400125628982]}, {\"name\": \"vT\", \"values\": [32.0069535572101]}, {\"name\": \"K_t\", \"values\": [5.943505949614364]}, {\"name\": \"K_i\", \"values\": [150.89852020829184]}, {\"name\": \"K_e\", \"values\": [5.90670286335029]}, {\"name\": \"Na_i\", \"values\": [12.702309626201476]}, {\"name\": \"Na_t\", \"values\": [132.39316813125978]}, {\"name\": \"Na_e\", \"values\": [132.99816117889287]}, {\"name\": \"n\", \"values\": [0.6110086409474482]}, {\"name\": \"h_K\", \"values\": [0.9679442466363977]}, {\"name\": \"m\", \"values\": [0.9989395956165837]}, {\"name\": \"h\", \"values\": [0.011749768142590286]}, {\"name\": \"S\", \"values\": [0.5803854365011085]}, {\"name\": \"n_t\", \"values\": [0.31830761662008844]}, {\"name\": \"h_K_t\", \"values\": [0.00497357422926022]}, {\"name\": \"m_t\", \"values\": [0.9991117400268303]}, {\"name\": \"h_t\", \"values\": [0.037558451115072077]}, {\"name\": \"S_t\", \"values\": [0.5807554861436337]}, {\"name\": \"O_0\", \"values\": [1.040185730088509e-06]}, {\"name\": \"O_1\", \"values\": [1.7764261016973102e-05]}, {\"name\": \"O_2\", \"values\": [0.00010214492966850287]}, {\"name\": \"O_3\", \"values\": [0.0001806534304265151]}, {\"name\": \"O_4\", \"values\": [7.488219517482619e-05]}, {\"name\": \"C_0\", \"values\": [0.5235895401090125]}, {\"name\": \"C_1\", \"values\": [0.3676938153420917]}, {\"name\": \"C_2\", \"values\": [0.0967148096280458]}, {\"name\": \"C_3\", \"values\": [0.011171890150204495]}, {\"name\": \"C_4\", \"values\": [0.0004534597686249599]}, {\"name\": \"Ca_1\", \"values\": [83.4425093078456]}, {\"name\": \"Ca_SR1\", \"values\": [1060.6659068598588]}, {\"name\": \"Ca_2\", \"values\": [3.0314298736642855]}, {\"name\": \"Ca_SR2\", \"values\": [1199.7566783570803]}, {\"name\": \"Ca_T_2\", \"values\": [24.733571546894552]}, {\"name\": \"Ca_P1\", \"values\": [641.6053654748897]}, {\"name\": \"Ca_P2\", \"values\": [616.2261077701777]}, {\"name\": \"Mg_P1\", \"values\": [810.8520611234248]}, {\"name\": \"Mg_P2\", \"values\": [810.9433332103665]}, {\"name\": \"Ca_Cs1\", \"values\": [16897.537754304743]}, {\"name\": \"Ca_Cs2\", \"values\": [16898.598304725354]}, {\"name\": \"Ca_ATP1\", \"values\": [108.46786460987256]}, {\"name\": \"Ca_ATP2\", \"values\": [10.511248908824985]}, {\"name\": \"Mg_ATP1\", \"values\": [7235.085496468541]}, {\"name\": \"Mg_ATP2\", \"values\": [7237.758581434536]}, {\"name\": \"ATP1\", \"values\": [656.44663892157]}, {\"name\": \"ATP2\", \"values\": [751.7301696566033]}, {\"name\": \"Mg1\", \"values\": [963.0444607968487]}, {\"name\": \"Mg2\", \"values\": [962.3184690076832]}, {\"name\": \"Ca_CaT2\", \"values\": [2.9818968299099726]}, {\"name\": \"D_0\", \"values\": [0.7991173574101228]}, {\"name\": \"D_1\", \"values\": [1.211475056337586]}, {\"name\": \"D_2\", \"values\": [2.9720138760087482]}, {\"name\": \"A_1\", \"values\": [0.29493676390460843]}, {\"name\": \"A_2\", \"values\": [0.23169305229112158]}, {\"name\": \"P\", \"values\": [0.2300207028869685]}, {\"name\": \"P_SR\", \"values\": [0.2301171644715283]}, {\"name\": \"P_C_SR\", \"values\": [0.22988283186639114]}, {\"name\": \"\", \"values\": [4.6429585851491e-310]}]}], \"timeStepNo\": 90001, \"currentTime\": 0.90001}";
  assertFileMatchesContent("out_0000009.py", referenceOutput);
}


TEST(CellMLTest, FastFibers)
{
  std::string pythonConfig = R"(

import numpy as np

# timing parameters
stimulation_frequency = 10.0      # [1/ms] frequency if which stimulation current can be switched on and off
dt_0D = 2e-4                      # timestep width of ODEs, cellml integration
dt_1D = 2e-3
dt_splitting = 2e-3
dt_3D = 0.5
output_timestep = 0.5
end_time = 30.0
n_elements = 100

stimulation_frequency = 100*1e-3   # [Hz]*1e-3 = [ms^-1]
frequency_jitter = 0.1
call_enable_begin = 15.0  # [s]*1e3 = [ms]

fiber_distribution_file = "../input/MU_fibre_distribution_10MUs.txt"
firing_times_file = "../input/MU_firing_times_always.txt"

# load MU distribution and firing times
fiber_distribution = np.genfromtxt(fiber_distribution_file, delimiter=" ")
firing_times = np.genfromtxt(firing_times_file)

def fiber_gets_stimulated(fiber_no, frequency, current_time):
  """
  determine if fiber fiber_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = 0

  # determine if fiber fires now
  index = int(np.round(current_time * frequency))
  n_firing_times = np.size(firing_times,0)

  #if firing_times[index % n_firing_times, mu_no] == 1:
  #print("fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(fiber_no, mu_no, current_time, (index % n_firing_times), firing_times[index % n_firing_times, mu_no], "true" if firing_times[index % n_firing_times, mu_no] == 1 else "false"))
  print("fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(fiber_no, mu_no, current_time, (index % n_firing_times), firing_times[index % n_firing_times, mu_no], "true" if firing_times[index % n_firing_times, mu_no] == 1 else "false"))

  return firing_times[index % n_firing_times, mu_no] == 1

# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, fiber_no):

  #print("call set_specific_states at time {}".format(current_time))

  # determine if fiber gets stimulated at the current time
  is_fiber_gets_stimulated = fiber_gets_stimulated(fiber_no, stimulation_frequency, current_time)

  if is_fiber_gets_stimulated:
    # determine nodes to stimulate (center node, left and right neighbour)
    #innervation_zone_width_n_nodes = innervation_zone_width*100  # 100 nodes per cm
    innervation_node_global = int(n_nodes_global / 2)  # + np.random.randint(-innervation_zone_width_n_nodes/2,innervation_zone_width_n_nodes/2+1)
    nodes_to_stimulate_global = [innervation_node_global]
    if innervation_node_global > 0:
      nodes_to_stimulate_global.insert(0, innervation_node_global-1)
    if innervation_node_global < n_nodes_global-1:
      nodes_to_stimulate_global.append(innervation_node_global+1)
    print("t: {}, stimulate fiber {} at nodes {}".format(current_time, fiber_no, nodes_to_stimulate_global))

    for node_no_global in nodes_to_stimulate_global:
      states[(node_no_global,0,0)] = 20.0   # key: ((x,y,z),nodal_dof_index,state_no)


# define the config dict
config = {
  "scenarioName": "not",
  "Meshes": {
    "MeshFiber_0": {
      "nElements": [n_elements],
      "physicalExtent": [n_elements/100.],
      "inputMeshIsGlobal": True,
    }
  },
  "Solvers": {
    "implicitSolver": {     # solver for the implicit timestepping scheme of the diffusion time step
      "maxIterations":      1e4,
      "relativeTolerance":  1e-10,
      "dumpFormat": "",
      "dumpFilename": "",
      "solverType": "gmres",
      "preconditionerType": "none"
    },
  },
  "RepeatedCall": {
    "timeStepWidth":          dt_splitting,
    "timeStepOutputInterval": 100,
    "endTime":                end_time,
    "MultipleInstances": {
      "ranksAllComputedInstances":  [0],
      "nInstances":                 1,
      "instances":
      [{
        "ranks": [0],
        "StrangSplitting": {
          #"numberTimeSteps": 1,
          "timeStepWidth":          dt_splitting,
          "timeStepOutputInterval": 100,
          "endTime":                dt_splitting,
          "transferSlotName":       "states",   # which output slot of the cellml adapter ("states" or "intermediates") to use for transfer to diffusion, in this case we need "states", states[0] which is Vm

          "Term1": {      # CellML, i.e. reaction term of Monodomain equation
            "MultipleInstances": {
              "logKey":             "duration_subdomains_z",
              "nInstances":         1,
              "instances":
              [{
                "ranks":                          [0],    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                "Heun" : {
                  "timeStepWidth":                dt_0D,  # 5e-5
                  "logTimeStepWidthAsKey":        "dt_0D",
                  "durationLogKey":               "duration_0D",
                  "initialValues":                [],
                  "timeStepOutputInterval":       1e4,
                  "inputMeshIsGlobal":            True,
                  "dirichletBoundaryConditions":  {},

                  "CellML" : {
                    "sourceFilename":                         "../input/hodgkin_huxley_1952.c",                          # input C++ source file, can be either generated by OpenCMISS or OpenCOR from cellml model
                    "compilerFlags":                          "-fPIC -O3 -shared ",
                    "useGivenLibrary":                        False,
                    "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                    "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                    "setSpecificStatesCallFrequency":         stimulation_frequency,   # set_specific_states should be called variables.stimulation_frequency times per ms
                    "setSpecificStatesFrequencyJitter":       frequency_jitter, # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                    "setSpecificStatesRepeatAfterFirstCall":  0.1,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                    "setSpecificStatesCallEnableBegin":       call_enable_begin,# [ms] first time when to call setSpecificStates
                    "additionalArgument":                     0,
                    "stimulationLogFilename":                 "out/stimulation_log.txt",

                    "outputIntermediateIndex":                0,                                              # which intermediate value to use in further computation
                    "outputStateIndex":                       0,                                              # Shorten / Hodgkin Huxley: state 0 = Vm, Shorten: rate 28 = gamma, intermediate 0 = gamma (OC_WANTED[0])
                    "parametersUsedAsIntermediate":           [],      #[32],       # list of intermediate value indices, that will be set by parameters. Explicitely defined parameters that will be copied to intermediates, this vector contains the indices of the algebraic array. This is ignored if the input is generated from OpenCMISS generated c code.
                    "parametersUsedAsConstant":               [2],          #[65],           # list of constant value indices, that will be set by parameters. This is ignored if the input is generated from OpenCMISS generated c code.
                    "parametersInitialValues":                [0.0],            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                    "meshName":                               "MeshFiber_0",
                    "prefactor":                              1.0,
                  },
                },
              }],
            }
          },
          "Term2": {     # Diffusion
            "MultipleInstances": {
              "nInstances": 1,
              "instances":
              [{
                "ranks":                         [0],    # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                "ImplicitEuler" : {
                  "initialValues":               [],
                  #"numberTimeSteps":            1,
                  "timeStepWidth":               dt_1D,  # 1e-5
                  "logTimeStepWidthAsKey":       "dt_1D",
                  "durationLogKey":              "duration_1D",
                  "timeStepOutputInterval":      1e4,
                  "dirichletBoundaryConditions": {}, #{0: -75.0036, -1: -75.0036},
                  "inputMeshIsGlobal":           True,
                  "solverName":                  "implicitSolver",
                  "FiniteElementMethod" : {
                    "maxIterations":             1e4,
                    "relativeTolerance":         1e-10,
                    "inputMeshIsGlobal":         True,
                    "meshName":                  "MeshFiber_0",
                    "prefactor":                 0.03,  # resolves to Conductivity / (Am * Cm)
                    "solverName":                "implicitSolver",
                  },
                  "OutputWriter" : [
                    #{"format": "Paraview", "outputInterval": int(1./variables.dt_1D*variables.output_timestep), "filename": "out/fiber_"+str(fiber_no), "binary": True, "fixedFormat": False, "combineFiles": True},
                    #{"format": "Paraview", "outputInterval": 1./variables.dt_1D*variables.output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                    #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "sphereSize": "0.02*0.02*0.02"},
                    #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "binary":True, "onlyNodalValues":True},
                  ]
                },
              }],
              "OutputWriter" : [
                {"format": "PythonFile", "outputInterval": int(1./dt_splitting*output_timestep), "filename": "out/not_fast/fibers", "binary": True, "fixedFormat": False, "combineFiles": True}
              ]
            },
          },
        }
      }]
    },
    "fiberDistributionFile":    fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
    "firingTimesFile":          firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
  }
}

)";

  DihuContext settings1(argc, argv, pythonConfig);


  // define problem
  TimeSteppingScheme::RepeatedCall<
    Control::MultipleInstances<                       // fibers
      OperatorSplitting::Strang<
        Control::MultipleInstances<
          TimeSteppingScheme::Heun<                   // fiber reaction term
            CellmlAdapter<
              4, 9,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
              FunctionSpace::FunctionSpace<
                Mesh::StructuredDeformableOfDimension<1>,
                BasisFunction::LagrangeOfOrder<1>
              >
            >
          >
        >,
        Control::MultipleInstances<
          TimeSteppingScheme::ImplicitEuler<          // fiber diffusion, note that implicit euler gives lower error in this case than crank nicolson
            SpatialDiscretization::FiniteElementMethod<
              Mesh::StructuredDeformableOfDimension<1>,
              BasisFunction::LagrangeOfOrder<1>,
              Quadrature::Gauss<2>,
              Equation::Dynamic::IsotropicDiffusion
            >
          >
        >
      >
    >
  > problem1(settings1);

  problem1.run();

  std::string strToReplace("out/not_fast/fibers");
  std::size_t pos = pythonConfig.find(strToReplace);
  pythonConfig.replace(pos, strToReplace.length(), "out/fast/fibers");

  DihuContext settings2(argc, argv, pythonConfig);

  // define problem with FastMonodomainSolver
  TimeSteppingScheme::RepeatedCall<
    FastMonodomainSolver<                        // a wrapper that improves performance of multidomain
      Control::MultipleInstances<                       // fibers
        OperatorSplitting::Strang<
          Control::MultipleInstances<
            TimeSteppingScheme::Heun<                   // fiber reaction term
              CellmlAdapter<
                4, 9,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
                FunctionSpace::FunctionSpace<
                  Mesh::StructuredDeformableOfDimension<1>,
                  BasisFunction::LagrangeOfOrder<1>
                >
              >
            >
          >,
          Control::MultipleInstances<
            TimeSteppingScheme::ImplicitEuler<          // fiber diffusion, note that implicit euler gives lower error in this case than crank nicolson
              SpatialDiscretization::FiniteElementMethod<
                Mesh::StructuredDeformableOfDimension<1>,
                BasisFunction::LagrangeOfOrder<1>,
                Quadrature::Gauss<2>,
                Equation::Dynamic::IsotropicDiffusion
              >
            >
          >
        >
      >
    >
  > problem2(settings2);
  
  // run problem
  problem2.run();

  // compare results of problem1 and problem2
  std::string command = R"(
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import py_reader
import numpy as np

print(os.getcwd());
directory1 = "out/fast"
directory2 = "out/not_fast"

# read all files in directories
files1 = []
for filename in os.listdir(directory1):
  if filename.endswith(".py"):
    files1.append(os.path.join(directory1, filename))
files1 = sorted(files1)

files2 = []
for filename in os.listdir(directory2):
  if filename.endswith(".py"):
    files2.append(os.path.join(directory2, filename))
files2 = sorted(files2)

print("files: ",files1,files2)

# load data
data1 = py_reader.load_data(files1)
data2 = py_reader.load_data(files2)

n_values = min(len(data1), len(data2))

if len(data1) != len(data2):
  print("Warning: Directory {} contains {} files, directory {} contains {} files.".format(directory1, len(data1), directory2, len(data2)))

component_name = "0"
total_error = 0
for i in range(n_values):
  values1 = py_reader.get_values(data1[i], "solution", component_name)
  values2 = py_reader.get_values(data2[i], "solution", component_name)

  error = np.linalg.norm(values1-values2) / np.size(values1);
  total_error += error

  print("file no. {}, error: {}".format(i, error))

total_error /= n_values
print("avg error: {}".format(total_error))

)";
  int returnValue = PyRun_SimpleString(command.c_str());
  LOG(DEBUG) << returnValue;
  ASSERT_EQ(returnValue, 0);
  PythonUtility::checkForError();

  // load main module and extract config
  PyObject *mainModule = PyImport_AddModule("__main__");
  PyObject *totalError = PyObject_GetAttrString(mainModule, "total_error");

  double error = PythonUtility::convertFromPython<double>::get(totalError);
  LOG(DEBUG) << "error: " << error;

  ASSERT_LE(error, 0.1);
}

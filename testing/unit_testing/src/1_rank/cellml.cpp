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


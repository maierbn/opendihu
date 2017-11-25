#pragma once

#include <vector>

#include "control/runnable.h"
#include "control/dihu_context.h"
#include "time_stepping_scheme/discretizable_in_time.h"

class CellmlAdapter : public DiscretizableInTime
{
public:
 
  ///! constructor
  CellmlAdapter(const DihuContext &context);
  
  ///! destructor
  ~CellmlAdapter();
  
  ///! load model
  void initialize();
  
  ///! evaluate rhs
  void evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime);
  
  ///! register a callback function setParameters that can set parameter values before each computation
  void registerSetParameters(void (*setParameters) (void *context, int nInstances, int timeStepNo, double currentTime, std::vector<double> &parameters));
  
  ///! directly call the python callback if it exists
  void callPythonSetParametersFunction(int nInstances, int timeStepNo, double currentTime, std::vector<double> &parameters);
  
  ///! directly call the python callback if it exists
  void callPythonHandleResultFunction(int nInstances, int timeStepNo, double currentTime, double *states, double *intermediates);
  
  ///! register a callbackfunction handleResult that gets called after each new values are available
  void registerHandleResult(void (*handleResult) (void *context, int nInstances, int timeStepNo, double currentTime, 
                                                    double *states, double intermediates[]));
  
  ///! set initial values as given in python config
  bool setInitialValues(Vec& initialValues);
  
  //! return a the mesh
  std::shared_ptr<Mesh::Mesh> mesh();
  
  //! return the number of states per instance
  int numberDegreesOfFreedomPerNode();
  
  //! get number of states, number of instances, number of intermediates and number of parameters
  void getNumbers(int &nStates, int &nInstances, int &nIntermediates, int &nParameters);
  
  //! return false because the object is independent of mesh type
  bool knowsMeshType();
  
private:
 
  //! initialize the rhs routine, either directly from a library or compile it
  void initializeRhsRoutine();  
  
  //! given a normal cellml source file for rhs routine create a second file for multiple instances. @return: if successful
  bool createSimdSourceFile(std::string &simdSourceFilename);
  
  //! scan the given cellml source file for initial values that are given by dummy assignments
  bool scanInitialValues(std::string sourceFilename, std::vector<double> &statesInitialValues);
  
  const DihuContext &context_;    ///< the context object containing everything to be stored
  
  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  
  void (*rhsRoutine_)(double, double *, double *, double *, double *);   ///< function pointer to a rhs function that is passed as dynamic library, computes rates and intermediate values from states. The parameters are: VOI, STATES, RATES, WANTED, KNOWN, (VOI: unclear, what it means)
  void (*rhsRoutineSimd_)(void *context, double *, double *, double *, double *);  ///< same functionality as rhsRoutine, however, can compute several instances of the problem in parallel. Data is assumed to contain values for a state contiguously, e.g. (state[1], state[1], state[1], state[2], state[2], state[2], ...). The first parameter is a this pointer
  void (*setParameters_) (void *context, int nInstances, int timeStepNo, double currentTime, std::vector<double> &parmeters);  ///< callback function that will be called before new states are computed. It can set new parameters ("known" variables) for the computation.
  void (*handleResult_) (void *context, int nInstances, int timeStepNo, double currentTime, double *states, double *intermediates);   ///< callback function that will be called after new states and intermediates were computed 
  
  std::shared_ptr<Mesh::Mesh> mesh_;    ///< a mesh, there are as many instances of the same CellML problem as there are nodes in the mesh
  
  int nStates_;           ///< number of states in one instance of the CellML problem
  int nInstances_;        ///< number of instances of the CellML problem, equals number of mesh nodes
  int nParameters_;       ///< number of parameters (=CellML name "known") in one instance of the CellML problem
  int nIntermediates_;    ///< number of intermediate values (=CellML name "wanted") in one instance of the CellML problem
  int setParametersCallInterval_;      ///< setParameters_ will be called every callInterval_ time steps
  int handleResultCallInterval_;      ///< handleResult will be called every callInterval_ time steps
  PyObject *pythonSetParametersFunction_;   ///< Python function handle that is called to set parameters to the CellML problem from the python config
  PyObject *pythonHandleResultFunction_;   ///< Python function handle that is called to process results from CellML problem from the python config
  
  std::string sourceFilename_;        ///< the filename of the source file that actually is used for rhs
  
  //std::vector<double> states_;    ///< vector of states, that are computed by rhsRoutine
  //std::vector<double> rates_;     ///< vector of rates, that are computed by rhsRoutine
  std::vector<double> parameters_; ///< vector of values that will be provided to CellML by the code, given by python config, CellML name: known
  std::vector<double> intermediates_;    ///< vector of intermediate values in DAE system. These can be computed directly from the actual states at any time. Gets computed by rhsRoutine from states, together with rates. OpenCMISS name is intermediate, CellML name: wanted
};
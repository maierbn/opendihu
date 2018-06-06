#pragma once

#include <vector>

#include "control/runnable.h"
#include "control/dihu_context.h"
#include "time_stepping_scheme/discretizable_in_time.h"
#include "output_writer/manager.h"
#include "basis_on_mesh/basis_on_mesh.h"
#include "basis_function/lagrange.h"
#include "cellml/01_rhs_routine_handler.h"

/** The is a class that contains cellml equations and can be used with a time stepping scheme.
 *  The nStates template parameter specifies the number of state variables that should be used with the integrator.
 *  It is necessary that this value is fixed at compile time because the timestepping scheme needs to know which field variable types is has to construct.
 *  This class can also be computed easily in multiple instances along the nodes of a mesh.
 */
template <int nStates>
class CallbackHandler :
  public RhsRoutineHandler<nStates>
{
public:

  ///! constructor
  CallbackHandler(DihuContext context);

  ///! destructor
  virtual ~CallbackHandler();
  
  ///! register a callback function setParameters that can set parameter values before each computation
  void registerSetParameters(void (*setParameters) (void *context, int nInstances, int timeStepNo, double currentTime,
                                                    std::vector<double> &parameters));

  ///! register a callbackfunction handleResult that gets called after each new values are available
  void registerHandleResult(void (*handleResult) (void *context, int nInstances, int timeStepNo, double currentTime,
                                                  double *states, double intermediates[]));
  
  ///! directly call the python callback if it exists
  void callPythonSetParametersFunction(int nInstances, int timeStepNo, double currentTime, std::vector<double> &parameters);

  ///! directly call the python callback if it exists
  void callPythonHandleResultFunction(int nInstances, int timeStepNo, double currentTime, double *states, double *intermediates);

protected:
 
  //! construct the python call back functions from config
  virtual void initializeCallbackFunctions();

  void (*setParameters_) (void *context, int nInstances, int timeStepNo, double currentTime, std::vector<double> &parmeters);  ///< callback function that will be called before new states are computed. It can set new parameters ("known" variables) for the computation.
  void (*handleResult_) (void *context, int nInstances, int timeStepNo, double currentTime, double *states, double *intermediates);   ///< callback function that will be called after new states and intermediates were computed

  int setParametersCallInterval_;      ///< setParameters_ will be called every callInterval_ time steps
  int handleResultCallInterval_;      ///< handleResult will be called every callInterval_ time steps
 
  PyObject *pythonSetParametersFunction_;   ///< Python function handle that is called to set parameters to the CellML problem from the python config
  PyObject *pythonHandleResultFunction_;   ///< Python function handle that is called to process results from CellML problem from the python config

  PyObject *pySetParametersFunctionAdditionalParameter_;  ///< an additional python object that will be passed as last argument to the setParameters callback function
  PyObject *pyHandleResultFunctionAdditionalParameter_;   ///< an additional python object that will be passed as last argument to the handleResult callback function

};

#include "cellml/02_callback_handler.tpp"
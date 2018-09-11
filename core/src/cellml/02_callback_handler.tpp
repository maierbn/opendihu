#include "cellml/02_callback_handler.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/mesh_manager.h"


template<int nStates, typename FunctionSpaceType>
CallbackHandler<nStates,FunctionSpaceType>::
CallbackHandler(DihuContext context) :
  RhsRoutineHandler<nStates,FunctionSpaceType>(context),
  DiscretizableInTime(SolutionVectorMapping()),
  setParameters_(NULL), handleResult_(NULL),
  pythonSetParametersFunction_(NULL), pythonHandleResultFunction_(NULL),
  pySetParametersFunctionAdditionalParameter_(NULL), pyHandleResultFunctionAdditionalParameter_(NULL)
{
}

template<int nStates, typename FunctionSpaceType>
CallbackHandler<nStates,FunctionSpaceType>::
~CallbackHandler()
{
  Py_CLEAR(pythonSetParametersFunction_);
  Py_CLEAR(pythonHandleResultFunction_);
}

template<int nStates, typename FunctionSpaceType>
void CallbackHandler<nStates,FunctionSpaceType>::
initializeCallbackFunctions()
{
  if (PythonUtility::hasKey(this->specificSettings_, "setParametersFunction"))
  {
    pythonSetParametersFunction_ = PythonUtility::getOptionFunction(this->specificSettings_, "setParametersFunction");
    setParametersCallInterval_ = PythonUtility::getOptionInt(this->specificSettings_, "setParametersCallInterval", 1, PythonUtility::Positive);
    setParameters_ = [](void *context, int nInstances, int timeStepNo, double currentTime, std::vector<double> &parameters)
    {
      CallbackHandler *cellmlAdapter = (CallbackHandler *)context;
      cellmlAdapter->callPythonSetParametersFunction(nInstances, timeStepNo, currentTime, parameters);
    };
    LOG(DEBUG) << "registered setParameters function";
  }

  if (PythonUtility::hasKey(this->specificSettings_, "handleResultFunction"))
  {
    pythonHandleResultFunction_ = PythonUtility::getOptionFunction(this->specificSettings_, "handleResultFunction");
    handleResultCallInterval_ = PythonUtility::getOptionInt(this->specificSettings_, "handleResultCallInterval", 1, PythonUtility::Positive);
    handleResult_ = [](void *context, int nInstances, int timeStepNo, double currentTime, double *states, double *intermediates)
    {
      CallbackHandler *cellmlAdapter = (CallbackHandler *)context;
      cellmlAdapter->callPythonHandleResultFunction(nInstances, timeStepNo, currentTime, states, intermediates);
    };
    LOG(DEBUG) << "registered handleResult function";
  }
  
  if (PythonUtility::hasKey(this->specificSettings_, "setParametersFunctionAdditionalParameter"))
  {
    pySetParametersFunctionAdditionalParameter_ = PythonUtility::getOptionPyObject(this->specificSettings_, "setParametersFunctionAdditionalParameter");
  }
  else 
  {
    pySetParametersFunctionAdditionalParameter_ = Py_None;
  }
  
  if (PythonUtility::hasKey(this->specificSettings_, "handleResultFunctionAdditionalParameter"))
  {
    pyHandleResultFunctionAdditionalParameter_ = PythonUtility::getOptionPyObject(this->specificSettings_, "handleResultFunctionAdditionalParameter");
  }
  else 
  {
    pyHandleResultFunctionAdditionalParameter_ = Py_None;
  }
}

template<int nStates, typename FunctionSpaceType>
void CallbackHandler<nStates,FunctionSpaceType>::
callPythonSetParametersFunction(int nInstances, int timeStepNo, double currentTime, std::vector< double >& parameters)
{
  if (pythonSetParametersFunction_ == NULL)
    return;
  
  // compose callback function
  PyObject *parametersList = PythonUtility::convertToPythonList(parameters);
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O)", nInstances, timeStepNo, currentTime, parametersList, pySetParametersFunctionAdditionalParameter_);
  PyObject *returnValue = PyObject_CallObject(pythonSetParametersFunction_, arglist);

  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();

  // copy new values in parametersList to parameters_ vector
  for (unsigned int i=0; i<parameters.size(); i++)
  {
    PyObject *item = PyList_GetItem(parametersList, (Py_ssize_t)i);
    parameters[i] = PythonUtility::convertFromPython<double>::get(item);
  }

  // decrement reference counters for python objects
  Py_CLEAR(parametersList);
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

template<int nStates, typename FunctionSpaceType>
void CallbackHandler<nStates,FunctionSpaceType>::
callPythonHandleResultFunction(int nInstances, int timeStepNo, double currentTime,
                               double *states, double *intermediates)
{
  if (pythonHandleResultFunction_ == NULL)
    return;

  // compose callback function
  LOG(DEBUG) << "callPythonHandleResultFunction: nInstances: " << this->nInstances_<< ", nStates: " << nStates << ", nIntermediates: " << this->nIntermediates_;
  PyObject *statesList = PythonUtility::convertToPythonList(nStates*this->nInstances_, states);
  PyObject *intermediatesList = PythonUtility::convertToPythonList(this->nIntermediates_*this->nInstances_, intermediates);
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O,O)", nInstances, timeStepNo, currentTime, statesList, intermediatesList, pyHandleResultFunctionAdditionalParameter_);
  PyObject *returnValue = PyObject_CallObject(pythonHandleResultFunction_, arglist);

  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();

  // decrement reference counters for python objects
  Py_CLEAR(statesList);
  Py_CLEAR(intermediatesList);
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

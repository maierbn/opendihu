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
  DiscretizableInTime(),
  setParameters_(NULL), handleResult_(NULL),
  pythonSetParametersFunction_(NULL), pythonSetSpecificParametersFunction_(NULL), pythonHandleResultFunction_(NULL),
  pySetParametersFunctionAdditionalParameter_(NULL), pyHandleResultFunctionAdditionalParameter_(NULL), pyGlobalNaturalDofsList_(NULL)
{
}

template<int nStates, typename FunctionSpaceType>
CallbackHandler<nStates,FunctionSpaceType>::
~CallbackHandler()
{
  Py_CLEAR(pythonSetSpecificParametersFunction_);
  Py_CLEAR(pythonSetParametersFunction_);
  Py_CLEAR(pythonHandleResultFunction_);
}

template<int nStates, typename FunctionSpaceType>
void CallbackHandler<nStates,FunctionSpaceType>::
initializeCallbackFunctions()
{
  if (this->specificSettings_.hasKey("setParametersFunction"))
  {
    pythonSetParametersFunction_ = this->specificSettings_.getOptionFunction("setParametersFunction");
    setParametersCallInterval_ = this->specificSettings_.getOptionInt("setParametersCallInterval", 1, PythonUtility::Positive);
    setParameters_ = [](void *context, int nInstances, int timeStepNo, double currentTime, std::vector<double> &parameters)
    {
      CallbackHandler *cellmlAdapter = (CallbackHandler *)context;
      cellmlAdapter->callPythonSetParametersFunction(nInstances, timeStepNo, currentTime, parameters);
    };
    LOG(DEBUG) << "registered setParameters function";
    LOG(WARNING) << "You specified the \"setParametersFunction\" callback which is slow, consider using \"setSpecificParametersFunction\" instead!";
  }
  if (this->specificSettings_.hasKey("setSpecificParametersFunction"))
  {
    pythonSetSpecificParametersFunction_ = this->specificSettings_.getOptionFunction("setSpecificParametersFunction");
    setSpecificParametersCallInterval_ = this->specificSettings_.getOptionInt("setSpecificParametersCallInterval", 1, PythonUtility::Positive);
    setSpecificParameters_ = [](void *context, int nInstances, int timeStepNo, double currentTime, std::map<global_no_t,double> &globalParameters)
    {
      CallbackHandler *cellmlAdapter = (CallbackHandler *)context;
      cellmlAdapter->callPythonSetSpecificParametersFunction(nInstances, timeStepNo, currentTime, globalParameters);
    };
    LOG(DEBUG) << "registered setSpecificParameters function";
  }

  if (this->specificSettings_.hasKey("handleResultFunction"))
  {
    pythonHandleResultFunction_ = this->specificSettings_.getOptionFunction("handleResultFunction");
    handleResultCallInterval_ = this->specificSettings_.getOptionInt("handleResultCallInterval", 1, PythonUtility::Positive);
    handleResult_ = [](void *context, int nInstances, int timeStepNo, double currentTime, double *states, double *intermediates)
    {
      CallbackHandler *cellmlAdapter = (CallbackHandler *)context;
      cellmlAdapter->callPythonHandleResultFunction(nInstances, timeStepNo, currentTime, states, intermediates);
    };
    LOG(DEBUG) << "registered handleResult function";
  }
  
  if (this->specificSettings_.hasKey("setParametersFunctionAdditionalParameter"))
  {
    pySetParametersFunctionAdditionalParameter_ = this->specificSettings_.getOptionPyObject("setParametersFunctionAdditionalParameter");
  }
  else 
  {
    pySetParametersFunctionAdditionalParameter_ = Py_None;
  }
  
  if (this->specificSettings_.hasKey("handleResultFunctionAdditionalParameter"))
  {
    pyHandleResultFunctionAdditionalParameter_ = this->specificSettings_.getOptionPyObject("handleResultFunctionAdditionalParameter");
  }
  else 
  {
    pyHandleResultFunctionAdditionalParameter_ = Py_None;
  }
}

template<int nStates, typename FunctionSpaceType>
void CallbackHandler<nStates,FunctionSpaceType>::
callPythonSetParametersFunction(int nInstances, int timeStepNo, double currentTime, std::vector<double> &parameters)
{
  if (pythonSetParametersFunction_ == NULL)
    return;

  // compose callback function
  PyObject *globalParametersDict = PyDict_New();
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O,O)", this->functionSpace_->meshPartitionBase()->nDofsGlobal(),
                                    timeStepNo, currentTime, globalParametersDict, pySetParametersFunctionAdditionalParameter_);
  PyObject *returnValue = PyObject_CallObject(pythonSetParametersFunction_, arglist);

  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();

  // copy new values in parametersList to parameters_ vector
  // loop over dict entries
/*
    bool first = true;
    while (PyDict_Next(object, &pos, &key, &value))
    {
      if (!first)
        line << ",";
      first = false;

      line << std::endl << std::string(indent+2, ' ');

      if (PyUnicode_Check(key))
      {
        std::string keyString = pyUnicodeToString(key);
        line << keyString<< ": ";
      }
      else if (PyLong_Check(key))
      {
        std::string keyString = std::to_string(PyLong_AsLong(key));
        line << keyString<< ": ";
      }
      else
      {
        line << "(key is of unknown type): ";
      }

      line << getString(value, indent+2, 0);
    }
    line << std::endl << std::string(indent, ' ') << "}";
  }
  else
  {
    line << "<unknown type>";
  }
*/

  for (unsigned int i=0; i<parameters.size(); i++)
  {
    PyObject *item = PyList_GetItem(parametersList, (Py_ssize_t)i);
    parameters[i] = PythonUtility::convertFromPython<double>::get(item);
  }

  // decrement reference counters for python objects
  Py_CLEAR(globalParametersDict);
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

template<int nStates, typename FunctionSpaceType>
void CallbackHandler<nStates,FunctionSpaceType>::
callPythonSetSpecificParametersFunction(int nInstances, int timeStepNo, double currentTime, std::map<global_no_t,double> &globalParameters)
{
  if (pythonSetSpecificParametersFunction_ == NULL)
    return;

  // create list of global dof nos if it does not already exist
  if (pyGlobalNaturalDofsList_ == nullptr)
  {
    std::vector<global_no_t> dofNosGlobalNatural;
    this->functionSpace_->meshPartitionBase()->getDofNosGlobalNatural(dofNosGlobalNatural);
    pyGlobalNaturalDofsList_ = PythonUtility::convertToPythonList(dofNosGlobalNatural);
  }

  // compose callback function
  PyObject *parametersList = PythonUtility::convertToPythonList(parameters);
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O,O)", this->functionSpace_->meshPartitionBase()->nDofsGlobal(), timeStepNo, currentTime, parametersList, pyGlobalNaturalDofsList_, pySetParametersFunctionAdditionalParameter_);
  PyObject *returnValue = PyObject_CallObject(pythonSetSpecificParametersFunction_, arglist);

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

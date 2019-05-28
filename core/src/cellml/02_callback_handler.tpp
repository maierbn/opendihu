#include "cellml/02_callback_handler.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/mesh_manager.h"

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
CallbackHandler(DihuContext context, bool noNewOutputWriter) :
  RhsRoutineHandler<nStates,nIntermediates_,FunctionSpaceType>(context, noNewOutputWriter),
  DiscretizableInTime(),
  setParameters_(NULL), setSpecificParameters_(NULL), setSpecificStates_(NULL), handleResult_(NULL),
  pythonSetParametersFunction_(NULL), pythonSetSpecificParametersFunction_(NULL), pythonSetSpecificStatesFunction_(NULL), pythonHandleResultFunction_(NULL),
  pySetFunctionAdditionalParameter_(NULL), pyHandleResultFunctionAdditionalParameter_(NULL), pyGlobalNaturalDofsList_(NULL)
{
}

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
CallbackHandler(DihuContext context) :
  RhsRoutineHandler<nStates,nIntermediates_,FunctionSpaceType>(context),
  DiscretizableInTime(),
  setParameters_(NULL), setSpecificParameters_(NULL), setSpecificStates_(NULL), handleResult_(NULL),
  pythonSetParametersFunction_(NULL), pythonSetSpecificParametersFunction_(NULL), pythonSetSpecificStatesFunction_(NULL), pythonHandleResultFunction_(NULL),
  pySetFunctionAdditionalParameter_(NULL), pyHandleResultFunctionAdditionalParameter_(NULL), pyGlobalNaturalDofsList_(NULL)
{
}

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
~CallbackHandler()
{
  clearPyObjects();
}

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
void CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
clearPyObjects()
{
  // Py_CLEAR has no effect if the variable is NULL
  Py_CLEAR(pythonSetParametersFunction_);
  Py_CLEAR(pythonSetSpecificParametersFunction_);
  Py_CLEAR(pythonSetSpecificStatesFunction_);
  Py_CLEAR(pythonHandleResultFunction_);

  Py_CLEAR(pySetFunctionAdditionalParameter_);
  Py_CLEAR(pyHandleResultFunctionAdditionalParameter_);
  Py_CLEAR(pyGlobalNaturalDofsList_);
}

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
void CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
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

    pySetFunctionAdditionalParameter_ = this->specificSettings_.getOptionPyObject("additionalArgument", Py_None);
    LOG(DEBUG) << "registered setParameters function, call interval: " << setParametersCallInterval_;
    LOG(WARNING) << "You specified the \"setParametersFunction\" callback which is slow, consider using \"setSpecificParametersFunction\" instead!";
  }

  if (this->specificSettings_.hasKey("setSpecificParametersFunction"))
  {
    pythonSetSpecificParametersFunction_ = this->specificSettings_.getOptionFunction("setSpecificParametersFunction");
    setSpecificParametersCallInterval_ = this->specificSettings_.getOptionInt("setSpecificParametersCallInterval", 1, PythonUtility::Positive);
    setSpecificParameters_ = [](void *context, int nInstances, int timeStepNo, double currentTime, std::vector<double> &localParameters)
    {
      CallbackHandler *cellmlAdapter = (CallbackHandler *)context;
      cellmlAdapter->callPythonSetSpecificParametersFunction(nInstances, timeStepNo, currentTime, localParameters);
    };

    pySetFunctionAdditionalParameter_ = this->specificSettings_.getOptionPyObject("additionalArgument", Py_None);
    LOG(DEBUG) << "registered setSpecificParameters function, call interval: " << setSpecificParametersCallInterval_;
  }

  if (this->specificSettings_.hasKey("setSpecificStatesFunction"))
  {
    pythonSetSpecificStatesFunction_ = this->specificSettings_.getOptionFunction("setSpecificStatesFunction");
    setSpecificStatesCallInterval_ = this->specificSettings_.getOptionInt("setSpecificStatesCallInterval", 1, PythonUtility::NonNegative);
    setSpecificStates_ = [](void *context, int nInstances, int timeStepNo, double currentTime, double *localStates)
    {
      CallbackHandler *cellmlAdapter = (CallbackHandler *)context;
      cellmlAdapter->callPythonSetSpecificStatesFunction(nInstances, timeStepNo, currentTime, localStates);
    };

    pySetFunctionAdditionalParameter_ = this->specificSettings_.getOptionPyObject("additionalArgument", Py_None);
    LOG(DEBUG) << "registered setSpecificStates function, call interval: " << setSpecificStatesCallInterval_;
  }

  if (this->specificSettings_.hasKey("handleResultFunction"))
  {
    pythonHandleResultFunction_ = this->specificSettings_.getOptionFunction("handleResultFunction");
    handleResultCallInterval_ = this->specificSettings_.getOptionInt("handleResultCallInterval", 1, PythonUtility::Positive);
    handleResult_ = [](void *context, int nInstances, int timeStepNo, double currentTime, double *localStates, double *intermediates)
    {
      CallbackHandler *cellmlAdapter = (CallbackHandler *)context;
      cellmlAdapter->callPythonHandleResultFunction(nInstances, timeStepNo, currentTime, localStates, intermediates);
    };

    pyHandleResultFunctionAdditionalParameter_ = this->specificSettings_.getOptionPyObject("handleResultFunctionAdditionalParameter", Py_None);
    LOG(DEBUG) << "registered handleResult function, call interval: " << handleResultCallInterval_;
  }
  
  if (this->specificSettings_.hasKey("setParametersFunctionAdditionalParameter"))
  {
    LOG(ERROR) << "Settings key \"setParametersFunctionAdditionalParameter\" has recently been changed to \"additionalArgument\".";
  }
}

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
void CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
callPythonSetParametersFunction(int nInstances, int timeStepNo, double currentTime, std::vector<double> &parameters)
{
  if (pythonSetParametersFunction_ == NULL)
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
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O,O)", this->functionSpace_->meshPartitionBase()->nDofsGlobal(),
                                    timeStepNo, currentTime, parametersList, pyGlobalNaturalDofsList_, pySetFunctionAdditionalParameter_);
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


template<int nStates, int nIntermediates_, typename FunctionSpaceType>
void CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
callPythonSetSpecificParametersFunction(int nInstances, int timeStepNo, double currentTime, std::vector<double> &localParameters)
{
  if (pythonSetSpecificParametersFunction_ == NULL)
    return;

  VLOG(1) << "callPythonSetSpecificParametersFunction timeStepNo=" << timeStepNo;

  // compose callback function
  PyObject *globalParametersDict = PyDict_New();
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O)", this->functionSpace_->meshPartitionBase()->nDofsGlobal(),
                                    timeStepNo, currentTime, globalParametersDict, pySetFunctionAdditionalParameter_);
  PyObject *returnValue = PyObject_CallObject(pythonSetSpecificParametersFunction_, arglist);

  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();

  const int D = FunctionSpaceType::dim();
  PyObject *keyPy, *valuePy;
  Py_ssize_t pos = 0;

  // iterate over key-value pairs in the python dict
  while (PyDict_Next(globalParametersDict, &pos, &keyPy, &valuePy))
  {
    // extract the key of dict that was given in the python config
    typedef std::tuple<std::array<global_no_t,D>,int,int> KeyType;
    KeyType key = PythonUtility::convertFromPython<KeyType>::get(keyPy);   // tuple( [x,y,z], nodalDofIndex, parameterNo )
    std::array<global_no_t,D> coordinatesGlobal = std::get<0>(key);
    int nodalDofIndex = std::get<1>(key);
    int parameterNo = std::get<2>(key);

    // extract the value
    double value = PythonUtility::convertFromPython<double>::get(valuePy);

    VLOG(1) << "coordinatesGlobal: " << coordinatesGlobal << ", nodalDofIndex: " << nodalDofIndex
      << ", parameterNo: " << parameterNo << ", value: " << value;

    // tranform global coordinates to local dof no, this can be outside of the local domain
    bool isOnLocalDomain;
    dof_no_t dofNoLocal = this->functionSpace_->getDofNoLocal(coordinatesGlobal, nodalDofIndex, isOnLocalDomain);
    dof_no_t nDofsLocalWithoutGhosts = this->functionSpace_->nDofsLocalWithoutGhosts();

    VLOG(1) << "dofNoLocal: " << dofNoLocal << " is OnLocalDomain: " << isOnLocalDomain;

    if (isOnLocalDomain)   // if the given parameter values is for a dof inside the current domain
    {
      // set first parameter value to given value
      localParameters[parameterNo*nDofsLocalWithoutGhosts + dofNoLocal] = value;
    }
    VLOG(1) << "localParameters: " << localParameters;
  }

  // decrement reference counters for python objects
  Py_CLEAR(globalParametersDict);
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
void CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
callPythonSetSpecificStatesFunction(int nInstances, int timeStepNo, double currentTime, double *localStates)
{
  if (pythonSetSpecificStatesFunction_ == NULL)
    return;

  VLOG(1) << "callPythonSetSpecificStatesFunction timeStepNo=" << timeStepNo;

  // compose callback function
  PyObject *globalStatesDict = PyDict_New();
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O)", this->functionSpace_->meshPartitionBase()->nDofsGlobal(),
                                    timeStepNo, currentTime, globalStatesDict, pySetFunctionAdditionalParameter_);
  PyObject *returnValue = PyObject_CallObject(pythonSetSpecificStatesFunction_, arglist);

  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();

  const int D = FunctionSpaceType::dim();
  PyObject *keyPy, *valuePy;
  Py_ssize_t pos = 0;

  // iterate over key-value pairs in the python dict
  while (PyDict_Next(globalStatesDict, &pos, &keyPy, &valuePy))
  {
    // extract the key of dict that was given in the python config
    typedef std::tuple<std::array<global_no_t,D>,int,int> KeyType;
    KeyType key = PythonUtility::convertFromPython<KeyType>::get(keyPy);   // tuple( [x,y,z], nodalDofIndex, stateNo )
    std::array<global_no_t,D> coordinatesGlobal = std::get<0>(key);
    int nodalDofIndex = std::get<1>(key);
    int stateNo = std::get<2>(key);

    // extract the value
    double value = PythonUtility::convertFromPython<double>::get(valuePy);

    VLOG(1) << "coordinatesGlobal: " << coordinatesGlobal << ", nodalDofIndex: " << nodalDofIndex
      << ", stateNo: " << stateNo << ", value: " << value;

    // tranform global coordinates to local dof no, this can be outside of the local domain
    bool isOnLocalDomain;
    dof_no_t dofNoLocal = this->functionSpace_->getDofNoLocal(coordinatesGlobal, nodalDofIndex, isOnLocalDomain);
    dof_no_t nDofsLocalWithoutGhosts = this->functionSpace_->nDofsLocalWithoutGhosts();

    VLOG(1) << "dofNoLocal: " << dofNoLocal << " is OnLocalDomain: " << isOnLocalDomain;

    if (isOnLocalDomain)   // if the given parameter values is for a dof inside the current domain
    {
      LOG(DEBUG) << "after python setSpecificStates callback: set dof " << dofNoLocal << ", state " << stateNo << " to value " << value;

      // set first parameter value to given value
      localStates[stateNo*nDofsLocalWithoutGhosts + dofNoLocal] = value;
    }
  }

  // decrement reference counters for python objects
  Py_CLEAR(globalStatesDict);
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
void CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
callPythonHandleResultFunction(int nInstances, int timeStepNo, double currentTime,
                               double *localStates, double *intermediates)
{
  if (pythonHandleResultFunction_ == NULL)
    return;

  // compose callback function
  LOG(DEBUG) << "callPythonHandleResultFunction: nInstances: " << this->nInstances_<< ", nStates: " << nStates
    << ", nIntermediates: " << this->nIntermediates();
  PyObject *statesList = PythonUtility::convertToPythonList(nStates*this->nInstances_, localStates);
  PyObject *intermediatesList = PythonUtility::convertToPythonList(nIntermediates_*this->nInstances_, intermediates);
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

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
double CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
lastCallSpecificStatesTime()
{
  return this->lastCallSpecificStatesTime_;
}

template<int nStates, int nIntermediates_, typename FunctionSpaceType>
void CallbackHandler<nStates,nIntermediates_,FunctionSpaceType>::
setLastCallSpecificStatesTime(double lastCallSpecificStatesTime)
{
  this->lastCallSpecificStatesTime_ = lastCallSpecificStatesTime;
  LOG(DEBUG) << "now set lastCallSpecificStatesTime_ to " << lastCallSpecificStatesTime_;
}


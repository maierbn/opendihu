#include "cellml/02_callback_handler.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/mesh_manager/mesh_manager.h"

template<int nStates, int nAlgebraics_, typename FunctionSpaceType>
CallbackHandler<nStates,nAlgebraics_,FunctionSpaceType>::
CallbackHandler(DihuContext context, bool initializeOutputWriter) :
  RhsRoutineHandler<nStates,nAlgebraics_,FunctionSpaceType>(context, initializeOutputWriter),
  DiscretizableInTime(),
  pythonSetSpecificParametersFunction_(NULL), pythonSetSpecificStatesFunction_(NULL), pythonHandleResultFunction_(NULL),
  pySetFunctionAdditionalParameter_(NULL), pyHandleResultFunctionAdditionalParameter_(NULL), pyGlobalNaturalDofsList_(NULL)
{
}

template<int nStates, int nAlgebraics_, typename FunctionSpaceType>
CallbackHandler<nStates,nAlgebraics_,FunctionSpaceType>::
CallbackHandler(DihuContext context) :
  RhsRoutineHandler<nStates,nAlgebraics_,FunctionSpaceType>(context),
  DiscretizableInTime(),
  pythonSetSpecificParametersFunction_(NULL), pythonSetSpecificStatesFunction_(NULL), pythonHandleResultFunction_(NULL),
  pySetFunctionAdditionalParameter_(NULL), pyHandleResultFunctionAdditionalParameter_(NULL), pyGlobalNaturalDofsList_(NULL)
{
}

template<int nStates, int nAlgebraics_, typename FunctionSpaceType>
CallbackHandler<nStates,nAlgebraics_,FunctionSpaceType>::
~CallbackHandler()
{
  clearPyObjects();
}

template<int nStates, int nAlgebraics_, typename FunctionSpaceType>
void CallbackHandler<nStates,nAlgebraics_,FunctionSpaceType>::
clearPyObjects()
{
  // Py_CLEAR has no effect if the variable is NULL
  Py_CLEAR(pythonSetSpecificParametersFunction_);
  Py_CLEAR(pythonSetSpecificStatesFunction_);
  Py_CLEAR(pythonHandleResultFunction_);

  Py_CLEAR(pySetFunctionAdditionalParameter_);
  Py_CLEAR(pyHandleResultFunctionAdditionalParameter_);
  Py_CLEAR(pyGlobalNaturalDofsList_);
}

template<int nStates, int nAlgebraics_, typename FunctionSpaceType>
void CallbackHandler<nStates,nAlgebraics_,FunctionSpaceType>::
initializeCallbackFunctions()
{
  if (this->specificSettings_.hasKey("setParametersFunction"))
  {
    LOG(ERROR) << "Option " << this->specificSettings_ << "[\"setParametersFunction\"] is no longer supported. Use \"setSpecificParametersFunction\" instead.";
  }

  // parse the setSpecificParametersFunction callback function
  if (this->specificSettings_.hasKey("setSpecificParametersFunction"))
  {
    pythonSetSpecificParametersFunction_ = this->specificSettings_.getOptionFunction("setSpecificParametersFunction");

    // if a callback function was specified, also parse the call interval
    if (pythonSetSpecificParametersFunction_)
    {
      setSpecificParametersCallInterval_ = this->specificSettings_.getOptionInt("setSpecificParametersCallInterval", 1, PythonUtility::Positive);
      pySetFunctionAdditionalParameter_ = this->specificSettings_.getOptionPyObject("additionalArgument", Py_None);

      LOG(DEBUG) << "registered setSpecificParameters function, call interval: " << setSpecificParametersCallInterval_;
    }
  }

  if (this->specificSettings_.hasKey("setSpecificStatesFunction"))
  {
    pythonSetSpecificStatesFunction_ = this->specificSettings_.getOptionFunction("setSpecificStatesFunction");

    // if a callback function was specified, also parse the call interval
    if (pythonSetSpecificStatesFunction_)
    {
      setSpecificStatesCallInterval_ = this->specificSettings_.getOptionInt("setSpecificStatesCallInterval", 1, PythonUtility::NonNegative);
      pySetFunctionAdditionalParameter_ = this->specificSettings_.getOptionPyObject("additionalArgument", Py_None);

      LOG(DEBUG) << "registered setSpecificStates function, call interval: " << setSpecificStatesCallInterval_;
    }
  }

  if (this->specificSettings_.hasKey("handleResultFunction"))
  {
    pythonHandleResultFunction_ = this->specificSettings_.getOptionFunction("handleResultFunction");

    // if a callback function was specified, also parse the call interval
    if (pythonHandleResultFunction_)
    {
      handleResultCallInterval_ = this->specificSettings_.getOptionInt("handleResultCallInterval", 1, PythonUtility::Positive);
      pyHandleResultFunctionAdditionalParameter_ = this->specificSettings_.getOptionPyObject("handleResultFunctionAdditionalParameter", Py_None);

      LOG(DEBUG) << "registered handleResult function, call interval: " << handleResultCallInterval_;
    }
  }
}

template<int nStates, int nAlgebraics_, typename FunctionSpaceType>
void CallbackHandler<nStates,nAlgebraics_,FunctionSpaceType>::
callPythonSetSpecificParametersFunction(int nInstances, int timeStepNo, double currentTime, double *localParameters, int nParameters)
{
  if (pythonSetSpecificParametersFunction_ == NULL)
    return;

  VLOG(1) << "callPythonSetSpecificParametersFunction timeStepNo=" << timeStepNo;

  // compose callback function
  PyObject *globalParametersDict = PyDict_New();
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O)", this->functionSpace_->meshPartition()->nDofsGlobal(),
                                    timeStepNo, currentTime, globalParametersDict, pySetFunctionAdditionalParameter_);
  PyObject *returnValue = PyObject_CallObject(pythonSetSpecificParametersFunction_, arglist);

  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PythonUtility::checkForError();

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

    // error checking on parameter
    if (parameterNo >= this->cellmlSourceCodeGenerator_.nParameters())
    {
      LOG(FATAL) << "In setSpecificParametersFunction: the parameters have an assignment "
          << "parameters[(coordinatesGlobal=" << coordinatesGlobal << ", nodalDofIndex=" << nodalDofIndex
          << ", parameterNo=" << parameterNo <<")] = " << value << ". But there are only " << this->cellmlSourceCodeGenerator_.nParameters()
          << " specified. Set \"parametersUsedAsAlgebraic\" and \"parametersUsedAsConstant\" appropriately.";
    }

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
      const int index = parameterNo*nDofsLocalWithoutGhosts + dofNoLocal;

      if (index >= nParameters*nDofsLocalWithoutGhosts)
      {
        LOG(FATAL) << "in setSpecificParametersFunction: the parameters have an assignment "
          << "parameters[(coordinatesGlobal=" << coordinatesGlobal << ", nodalDofIndex=" << nodalDofIndex
          << ", parameterNo=" << parameterNo <<")] = " << value << ".\n "
          << "The global coordinates, " << coordinatesGlobal << " refer to local dof " << dofNoLocal
          << ". There are " << nDofsLocalWithoutGhosts << " local dofs in total.\n"
          << "There are " << this->cellmlSourceCodeGenerator_.nParameters() << "=" << nParameters << " parameters.\n"
          << "The local field of parameters has " << nParameters*nDofsLocalWithoutGhosts << " entries, however the index is "
          << parameterNo << "*" << nDofsLocalWithoutGhosts << " + " << dofNoLocal << " = " << index
          << " >= " << nParameters*nDofsLocalWithoutGhosts;
      }

      // the local parameter array stores the parameter values in struct-of-array layout, there is space allocated for `nAlgebraics` parameters, even if there are only nLocalParameters parameters
      *(localParameters + index) = value;
    }
  }

  // decrement reference counters for python objects
  Py_CLEAR(globalParametersDict);
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

template<int nStates, int nAlgebraics_, typename FunctionSpaceType>
void CallbackHandler<nStates,nAlgebraics_,FunctionSpaceType>::
callPythonSetSpecificStatesFunction(int nInstances, int timeStepNo, double currentTime, double *localStates)
{
  if (pythonSetSpecificStatesFunction_ == NULL)
    return;

  VLOG(1) << "callPythonSetSpecificStatesFunction timeStepNo=" << timeStepNo;

  // compose callback function
  PyObject *globalStatesDict = PyDict_New();
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O)", this->functionSpace_->meshPartition()->nDofsGlobal(),
                                    timeStepNo, currentTime, globalStatesDict, pySetFunctionAdditionalParameter_);
  PyObject *returnValue = PyObject_CallObject(pythonSetSpecificStatesFunction_, arglist);

  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PythonUtility::checkForError();

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

template<int nStates, int nAlgebraics_, typename FunctionSpaceType>
void CallbackHandler<nStates,nAlgebraics_,FunctionSpaceType>::
callPythonHandleResultFunction(int nInstances, int timeStepNo, double currentTime,
                               double *localStates, double *algebraics)
{
  if (pythonHandleResultFunction_ == NULL)
    return;

  // compose callback function
  LOG(DEBUG) << "callPythonHandleResultFunction: nInstances: " << this->nInstances_<< ", nStates: " << nStates
    << ", nAlgebraics: " << this->nAlgebraics();
  PyObject *statesList = PythonUtility::convertToPythonList(nStates*this->nInstances_, localStates);
  PyObject *algebraicsList = PythonUtility::convertToPythonList(nAlgebraics_*this->nInstances_, algebraics);
  PyObject *arglist = Py_BuildValue("(i,i,d,O,O,O)", nInstances, timeStepNo, currentTime, statesList, algebraicsList, pyHandleResultFunctionAdditionalParameter_);
  PyObject *returnValue = PyObject_CallObject(pythonHandleResultFunction_, arglist);

  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PythonUtility::checkForError();

  // decrement reference counters for python objects
  Py_CLEAR(statesList);
  Py_CLEAR(algebraicsList);
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

template<int nStates, int nAlgebraics_, typename FunctionSpaceType>
double CallbackHandler<nStates,nAlgebraics_,FunctionSpaceType>::
lastCallSpecificStatesTime()
{
  return this->lastCallSpecificStatesTime_;
}

template<int nStates, int nAlgebraics_, typename FunctionSpaceType>
void CallbackHandler<nStates,nAlgebraics_,FunctionSpaceType>::
setLastCallSpecificStatesTime(double lastCallSpecificStatesTime)
{
  this->lastCallSpecificStatesTime_ = lastCallSpecificStatesTime;
  LOG(DEBUG) << "now set lastCallSpecificStatesTime_ to " << lastCallSpecificStatesTime_;
}


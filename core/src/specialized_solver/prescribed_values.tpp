#include "specialized_solver/prescribed_values.h"

#include <omp.h>
#include <sstream>

#include "control/diagnostic_tool/solver_structure_visualizer.h"

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
PrescribedValues(DihuContext context) :
  Runnable(),
  ::TimeSteppingScheme::TimeSteppingScheme(context["PrescribedValues"]),
  data_(this->context_), pyCallbackAdditionalParameter_(Py_None), pyGlobalNaturalDofsList_(nullptr)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  // the additional parameter for callback functions
  pyCallbackAdditionalParameter_ = this->specificSettings_.getOptionPyObject("additionalArgument", Py_None);
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute time span of this method
  double timeSpan = this->endTime_ - this->startTime_;

  // output for debugging
  LOG(DEBUG) << "PrescribedValues::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    // in defined intervals (settings "timeStepOutputInterval") print out the current timestep
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "PrescribedValues, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // call the callback functions to update the field variable values
    callCallbacks(timeStepNo, currentTime);

    // advance simulation time
    timeStepNo++;

    // compute new current simulation time
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values using the output writers
    this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
initialize()
{
  // initialize() will be called before the simulation starts.

  // call initialize of the parent class, this parses the timestepping settings from the settings file
  TimeSteppingScheme::TimeSteppingScheme::initialize();

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("PrescribedValues", false);   // hasInternalConnectionToFirstNestedSolver=false (the last argument) means slot connector data is not shared with the first subsolver
  // if you have your own slot connector data rather than the one of the subsolver, call "addSolver" with false as second argument


  // retrieve the function space
  std::shared_ptr<FunctionSpaceType> functionSpace = this->context_.meshManager()->functionSpace<FunctionSpaceType>(specificSettings_);

  // Pass the function space to the data object. data_ stores field variables.
  // It needs to know the number of nodes and degrees of freedom (dof) from the function space in order to create the vectors with the right size.
  data_.setFunctionSpace(functionSpace);

  /** example for the settings of a field variable:
   *
   * fieldVariables1: [
   *  {"name": "a", "callback": callback},
   * ]
   */

  // parse field variables

  std::vector<std::string> fieldVariable1Names;
  std::vector<std::string> fieldVariable2Names;

  // parse option "fieldVariables1"
  std::vector<PyObject *> fieldVariable1Settings;
  this->specificSettings_.getOptionVector<PyObject *>("fieldVariables1", fieldVariable1Settings);

  // loop over field variable settings
  for (PyObject *pyObject : fieldVariable1Settings)
  {
    PythonConfig settings(this->specificSettings_, "fieldVariables1", pyObject);
    std::string name = settings.getOptionString("name", "a");
    PyObject *pyCallback = settings.getOptionPyObject("callback");

    if (pyCallback == Py_None)
    {
      fieldVariable1Names.push_back(name);
      callbackFunctions1_.push_back(Py_None);
    }
    else if (!PyFunction_Check(pyCallback))
    {
      LOG(FATAL) << "Given value under " << settings << "[\"callback\"] is not a python function.";
    }
    else
    {
      fieldVariable1Names.push_back(name);
      callbackFunctions1_.push_back(pyCallback);
    }
  }

  // parse option "fieldVariables2"
  std::vector<PyObject *> fieldVariable2Settings;
  this->specificSettings_.getOptionVector<PyObject *>("fieldVariables2", fieldVariable2Settings);

  // loop over field variable settings
  for (PyObject *pyObject : fieldVariable2Settings)
  {
    PythonConfig settings(this->specificSettings_, "fieldVariables2", pyObject);
    std::string name = settings.getOptionString("name", "b");
    PyObject *pyCallback = settings.getOptionPyObject("callback");

    if (pyCallback == Py_None)
    {
      fieldVariable2Names.push_back(name);
      callbackFunctions2_.push_back(Py_None);
    }
    else if (!PyFunction_Check(pyCallback))
    {
      LOG(FATAL) << "Given value under " << settings << "[\"callback\"] is not a python function.";
    }
    else
    {
      fieldVariable2Names.push_back(name);
      callbackFunctions2_.push_back(pyCallback);
    }
  }

  // now call initialize, data will then create all variables (Petsc Vec's)
  data_.initialize(fieldVariable1Names, fieldVariable2Names);

  // set the slotConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setSlotConnectorData(getSlotConnectorData());
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
callCallbacks(int timeStepNo, double currentTime)
{
  // call callback functions for field variable 1
  for (int fieldVariable1No = 0; fieldVariable1No < callbackFunctions1_.size(); fieldVariable1No++)
  {
    // do not call callback if the function is None
    if (callbackFunctions1_[fieldVariable1No] == Py_None)
      continue;

    assert (callbackFunctions1_[fieldVariable1No]);
    assert (PyFunction_Check(callbackFunctions1_[fieldVariable1No]));

    //callbackFunctions1_
    // get all values
    std::vector<double> values;
    this->data_.fieldVariable1(fieldVariable1No)->getValuesWithoutGhosts(values);

    // create list of global dof nos if it does not already exist
    if (pyGlobalNaturalDofsList_ == nullptr)
    {
      //! get a vector of global natural dof nos of the locally stored non-ghost dofs, needed for setParameters callback function in cellml adapter
      std::vector<global_no_t> dofNosGlobalNatural;
      this->data_.functionSpace()->meshPartition()->getDofNosGlobalNatural(dofNosGlobalNatural);
      pyGlobalNaturalDofsList_ = PythonUtility::convertToPythonList(dofNosGlobalNatural);
    }

    const int D = FunctionSpaceType::dim();
    std::array<long,D> nNodesGlobalPerCoordinateDirection;
    for (int coordinateDirection = 0; coordinateDirection < D; coordinateDirection++)
    {
      nNodesGlobalPerCoordinateDirection[coordinateDirection] = this->data_.functionSpace()->meshPartition()->nNodesGlobal(coordinateDirection);
    }

    PyObject *nNodesGlobalPerCoordinateDirectionPy = PythonUtility::convertToPythonList<D>(nNodesGlobalPerCoordinateDirection);

    // compose arguments for callback function
    PyObject *valuesList = PythonUtility::convertToPythonList(values);
    PyObject *arglist = Py_BuildValue("(i,O,i,d,O,O,O)", this->data_.functionSpace()->meshPartition()->nDofsGlobal(), nNodesGlobalPerCoordinateDirectionPy,
                                      timeStepNo, currentTime, valuesList, pyGlobalNaturalDofsList_, pyCallbackAdditionalParameter_);
    PyObject *returnValue = PyObject_CallObject(callbackFunctions1_[fieldVariable1No], arglist);

    // if there was an error while executing the function, print the error message
    if (returnValue == NULL)
      PyErr_Print();

    // copy new values in valuesList to field variable
    for (unsigned int i=0; i<values.size(); i++)
    {
      PyObject *item = PyList_GetItem(valuesList, (Py_ssize_t)i);
      values[i] = PythonUtility::convertFromPython<double>::get(item);
    }

    VLOG(1) << "set values to field variable1(index " << fieldVariable1No << ") \"" << this->data_.fieldVariable1(fieldVariable1No)->name()
      << ": " << values;

    // set the values in field variable
    this->data_.fieldVariable1(fieldVariable1No)->setValuesWithoutGhosts(values);

    // decrement reference counters for python objects
    Py_CLEAR(valuesList);
    Py_CLEAR(returnValue);
    Py_CLEAR(arglist);
  }


  // call callback functions for field variable 2
  for (int fieldVariable2No = 0; fieldVariable2No < callbackFunctions2_.size(); fieldVariable2No++)
  {
    // do not call callback if the function is None
    if (callbackFunctions2_[fieldVariable2No] == Py_None)
      continue;

    //callbackFunctions2_
    // get all values
    std::vector<double> values;
    this->data_.fieldVariable2(fieldVariable2No)->getValuesWithoutGhosts(values);

    // create list of global dof nos if it does not already exist
    if (pyGlobalNaturalDofsList_ == nullptr)
    {
      //! get a vector of global natural dof nos of the locally stored non-ghost dofs, needed for setParameters callback function in cellml adapter
      std::vector<global_no_t> dofNosGlobalNatural;
      this->data_.functionSpace()->meshPartition()->getDofNosGlobalNatural(dofNosGlobalNatural);
      pyGlobalNaturalDofsList_ = PythonUtility::convertToPythonList(dofNosGlobalNatural);
    }

    // compose arguments for callback function
    PyObject *valuesList = PythonUtility::convertToPythonList(values);
    PyObject *arglist = Py_BuildValue("(i,i,d,O,O,O)", this->data_.functionSpace()->meshPartition()->nDofsGlobal(),
                                      timeStepNo, currentTime, valuesList, pyGlobalNaturalDofsList_, pyCallbackAdditionalParameter_);
    PyObject *returnValue = PyObject_CallObject(callbackFunctions2_[fieldVariable2No], arglist);

    // if there was an error while executing the function, print the error message
    if (returnValue == NULL)
      PyErr_Print();

    // copy new values in valuesList to field variable
    for (unsigned int i=0; i<values.size(); i++)
    {
      PyObject *item = PyList_GetItem(valuesList, (Py_ssize_t)i);
      values[i] = PythonUtility::convertFromPython<double>::get(item);
    }

    // set the values in field variable
    this->data_.fieldVariable2(fieldVariable2No)->setValuesWithoutGhosts(values);

    // decrement reference counters for python objects
    Py_CLEAR(valuesList);
    Py_CLEAR(returnValue);
    Py_CLEAR(arglist);
  }
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
run()
{
  // The run method should not be changed. It is the method that gets called directly from the example main file.
  // If this solver itself is nested in other solvers or coupling schemes,
  // run() will not be called, but the enclosing solver will call initialize() and advanceTimeSpan().
  initialize();

  advanceTimeSpan();
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
reset()
{
  // "uninitialize" everything
  callbackFunctions1_.clear();
  callbackFunctions2_.clear();
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
typename PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::Data &PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
data()
{
  // get a reference to the data object
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the slot_connector_data_transfer class
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::shared_ptr<typename PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::SlotConnectorDataType> PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
getSlotConnectorData()
{
  return data_.getSlotConnectorData();
}

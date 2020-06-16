#include "specialized_solver/my_new_solver/my_new_timestepping_solver.h"

#include <omp.h>
#include <sstream>

template<typename TimeStepping>
MyNewTimesteppingSolver<TimeStepping>::
MyNewTimesteppingSolver(DihuContext context) :
  Runnable(),
  ::TimeSteppingScheme::TimeSteppingScheme(context["MyNewTimesteppingSolver"]),   // replace "MyNewTimesteppingSolver" by the name of your solver, this will be the key for the dict in settings
  timeSteppingScheme_(this->context_),
  data_(this->context_)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  // parse options
  int myOption = this->specificSettings_.getOptionInt("myOption", 1, PythonUtility::Positive);

  LOG(DEBUG) << "myOption: " << myOption;
}

template<typename TimeStepping>
void MyNewTimesteppingSolver<TimeStepping>::
advanceTimeSpan()
{
  // This method computes some time steps of the simulation by running a for loop over the time steps.
  // The number of steps, timestep width and current time are all set by the parent class, TimeSteppingScheme.
  // You shouldn't change too much in this method.

  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute time span of this method
  double timeSpan = this->endTime_ - this->startTime_;

  // output for debugging
  LOG(DEBUG) << "MyNewTimesteppingSolver::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    // in defined intervals (settings "timeStepOutputInterval") print out the current timestep
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "MyNewTimesteppingSolver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // Now, define what the solver does in this time step. Typically you want to execute the nested timestepping object.

    // Set timespan for timeSteppingScheme_, the nested solver should advance its simulation by our timeStepWidth_.
    // This, in turn, may lead to multiple timesteps in the timeSteppingScheme_.
    this->timeSteppingScheme_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    // advance the simulation by the specified time span
    timeSteppingScheme_.advanceTimeSpan();

    // probably do something more here, maybe in a separate method:
    //executeMyHelperMethod();

    // ...


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

template<typename TimeStepping>
void MyNewTimesteppingSolver<TimeStepping>::
initialize()
{
  // initialize() will be called before the simulation starts.

  // call initialize of the parent class, this parses the timestepping settings from the settings file
  TimeSteppingScheme::TimeSteppingScheme::initialize();

  // add this solver to the solvers diagram, which is a SVG file that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("MyNewTimesteppingSolver", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means output connector data is shared with the first subsolver
  // if you have your own output connector data rather than the one of the subsolver, call "addSolver" with false as second argument

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested timestepping solver
  timeSteppingScheme_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // In order to initialize the data object and actuall create all variables, we first need to assign a function space to the data object.
  // A function space object of type FunctionSpace<MeshType,BasisFunctionType> (see "function_space/function_space.h")
  // is an object that stores the mesh (e.g., all nodes and elements) as well as the basis function (e.g. linear Lagrange basis functions).
  // The timeSteppingScheme_ solver already created a function space that we should use. We already have a typedef "FunctionSpace" that is the class of timeSteppingScheme_'s function space type.
  std::shared_ptr<FunctionSpace> functionSpace = timeSteppingScheme_.data().functionSpace();

  // Pass the function space to the data object. data_ stores field variables.
  // It needs to know the number of nodes and degrees of freedom (dof) from the function space in order to create the vectors with the right size.
  data_.setFunctionSpace(functionSpace);

  // now call initialize, data will then create all variables (Petsc Vec's)
  data_.initialize();

  // set the outputConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setOutputConnectorData(getOutputConnectorData());

  // here is the space to initialize anything else that is needed for your solver
}

template<typename TimeStepping>
void MyNewTimesteppingSolver<TimeStepping>::
run()
{
  // The run method should not be changed. It is the method that gets called directly from the example main file.
  // If this solver itself is nested in other solvers or coupling schemes,
  // run() will not be called, but the enclosing solver will call initialize() and advanceTimeSpan().
  initialize();

  advanceTimeSpan();
}

template<typename TimeStepping>
void MyNewTimesteppingSolver<TimeStepping>::
reset()
{
  timeSteppingScheme_.reset();

  // "uninitialize" everything
}

template<typename TimeStepping>
void MyNewTimesteppingSolver<TimeStepping>::
executeMyHelperMethod()
{
  // this is the template for an own private method

  // for example you can get the solution values of the timeSteppingScheme_ by
  Vec solution = timeSteppingScheme_.data().solution()->valuesGlobal();

  // As an example we invert all the solution values.
  // Because "Vec"'s are actually pointers, this effects the actual data, no copy-back is needed.
  // Note the error handling with Petsc functions. Always use "ierr = PetscFunction(); CHKERRV(ierr);"
  PetscErrorCode ierr;
  ierr = VecScale(solution, -1); CHKERRV(ierr);

  // get a field variable from data object:
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace,1>> fieldVariableA = data_.fieldVariableA();

  // copy the solution from the timeSteppingScheme_ to fieldVariableA.
  // note, you get the Petsc "Vec" of a field variable by "valuesGlobal()"
  ierr = VecCopy(solution, fieldVariableA->valuesGlobal()); CHKERRV(ierr);

  // add 1.0 to every entry in fieldVariableA, this already updates fieldVariableA in data because it is a pointer. There is no need to copy the values back.
  ierr = VecShift(fieldVariableA->valuesGlobal(), 1.0); CHKERRV(ierr);
}

template<typename TimeStepping>
typename MyNewTimesteppingSolver<TimeStepping>::Data &MyNewTimesteppingSolver<TimeStepping>::
data()
{
  // get a reference to the data object
  return data_;

  // The timeSteppingScheme_ object also has a data object, we could also directly use this and avoid having an own data object:
  //  return timeSteppingScheme_.data();
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class
template<typename TimeStepping>
std::shared_ptr<typename MyNewTimesteppingSolver<TimeStepping>::OutputConnectorDataType> MyNewTimesteppingSolver<TimeStepping>::
getOutputConnectorData()
{
  //! This is relevant only, if this solver is part of a splitting or coupling scheme. Then this method returns the values/variables that will be
  // transferred to the other solvers. We can just reuse the values of the timeSteppingScheme_.
  return timeSteppingScheme_.getOutputConnectorData();
}

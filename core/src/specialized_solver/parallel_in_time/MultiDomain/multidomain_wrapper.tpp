#include "specialized_solver/parallel_in_time/MultiDomain/multidomain_wrapper.h"

#include <omp.h>
#include <sstream>

template<typename StrangSplittingMultidomain>
MultidomainWrapper<StrangSplittingMultidomain>::
MultidomainWrapper(DihuContext context) :
  Runnable(),
  strangSplittingMultidomain_(context),
  context_(context), specificSettings_(context_.getPythonConfig())
{
}

template<typename StrangSplittingMultidomain>
void MultidomainWrapper<StrangSplittingMultidomain>::
advanceTimeSpan()
{
  strangSplittingMultidomain_.advanceTimeSpan();
}

template<typename StrangSplittingMultidomain>
void MultidomainWrapper<StrangSplittingMultidomain>::
initialize()
{
  outputWriterManager_.initialize(this->context_, this->specificSettings_);

  strangSplittingMultidomain_.initialize();
}

template<typename StrangSplittingMultidomain>
void MultidomainWrapper<StrangSplittingMultidomain>::
run()
{
  strangSplittingMultidomain_.run();
}

//! setup a new system matrix within the multidomain solver
template<typename StrangSplittingMultidomain>
void MultidomainWrapper<StrangSplittingMultidomain>::
setSystemMatrix(double timeStepWidth)
{
  // get the multidomain solver
  typename StrangSplittingMultidomain::TimeStepping2Type &multidomainSolver = strangSplittingMultidomain_.timeStepping2();

  // assemble a new system matrix with the given time step width, false = disable warning
  multidomainSolver.updateSystemMatrix(timeStepWidth, false);
}

//! get the local number of solution values, which is the number of entries set by getSolution
template<typename StrangSplittingMultidomain>
int MultidomainWrapper<StrangSplittingMultidomain>::
nSolutionValuesLocal()
{
  using HeunType = typename StrangSplittingMultidomain::TimeStepping1Type::TimeSteppingSchemeType;
  using CellMLAdapterType = typename HeunType::DiscretizableInTime;

  // get the first operator of the splitting
  typename StrangSplittingMultidomain::TimeStepping1Type &multipleInstancesHeun = strangSplittingMultidomain_.timeStepping1();
  assert(!multipleInstancesHeun.instancesLocal().empty());

  int nStates = CellMLAdapterType::nStates();
  int nCompartments = multipleInstancesHeun.instancesLocal().size();

  HeunType &heunScheme = multipleInstancesHeun.instancesLocal().front();
  CellMLAdapterType &cellmlAdapter = heunScheme.discretizableInTime();

  int nDofsLocalWithoutGhosts = cellmlAdapter.data().functionSpace()->nDofsLocalWithoutGhosts();
  return nDofsLocalWithoutGhosts * nStates * nCompartments;
}

template<typename StrangSplittingMultidomain>
void MultidomainWrapper<StrangSplittingMultidomain>::
getSolution(double *result, int timestepNo, double currentTime)
{
  // get the first operator of the splitting
  typename StrangSplittingMultidomain::TimeStepping1Type &multipleInstancesHeun = strangSplittingMultidomain_.timeStepping1();

  using HeunType = typename StrangSplittingMultidomain::TimeStepping1Type::TimeSteppingSchemeType;
  using CellMLAdapterType = typename HeunType::DiscretizableInTime;

  int nStates = CellMLAdapterType::nStates();
  int nCompartments = multipleInstancesHeun.instancesLocal().size();

  // determine number of local dofs
  assert(!multipleInstancesHeun.instancesLocal().empty());
  HeunType &heunScheme = multipleInstancesHeun.instancesLocal().front();
  CellMLAdapterType &cellmlAdapter = heunScheme.discretizableInTime();

  int nDofsLocalWithoutGhosts = cellmlAdapter.data().functionSpace()->nDofsLocalWithoutGhosts();

  std::vector<double> localValues;

  // loop over multidomain compartments
  for (int compartmentNo = 0; compartmentNo < nCompartments; compartmentNo++)
  {
    HeunType &heunScheme = multipleInstancesHeun.instancesLocal()[compartmentNo];
    CellMLAdapterType &cellmlAdapter = heunScheme.discretizableInTime();

    // loop over CellML states
    for (int stateNo = 0; stateNo < nStates; stateNo++)
    {
      // get the field variable entries for the current state of the current compartment
      localValues.clear();
      cellmlAdapter.data().states()->getValuesWithoutGhosts(stateNo, localValues);

      assert(nDofsLocalWithoutGhosts == localValues.size());

      // store the values in the result vector
      for (dof_no_t dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts; dofNoLocal++)
      {
        result[dofNoLocal * nCompartments * nStates + compartmentNo * nStates + stateNo] = localValues[dofNoLocal];
      }
    }
  }

  // call output writer to write output files
  this->outputWriterManager_.writeOutput(cellmlAdapter.data(), timestepNo, currentTime);
}

template<typename StrangSplittingMultidomain>
void MultidomainWrapper<StrangSplittingMultidomain>::
setSolution(double *data)
{
  // get the first operator of the splitting
  typename StrangSplittingMultidomain::TimeStepping1Type &multipleInstancesHeun = strangSplittingMultidomain_.timeStepping1();

  using HeunType = typename StrangSplittingMultidomain::TimeStepping1Type::TimeSteppingSchemeType;
  using CellMLAdapterType = typename HeunType::DiscretizableInTime;

  int nStates = CellMLAdapterType::nStates();
  int nCompartments = multipleInstancesHeun.instancesLocal().size();

  // determine number of local dofs
  assert(!multipleInstancesHeun.instancesLocal().empty());
  HeunType &heunScheme = multipleInstancesHeun.instancesLocal().front();
  CellMLAdapterType &cellmlAdapter = heunScheme.discretizableInTime();

  int nDofsLocalWithoutGhosts = cellmlAdapter.data().functionSpace()->nDofsLocalWithoutGhosts();

  std::vector<double> localValues(nDofsLocalWithoutGhosts);

  // loop over multidomain compartments
  for (int compartmentNo = 0; compartmentNo < nCompartments; compartmentNo++)
  {
    HeunType &heunScheme = multipleInstancesHeun.instancesLocal()[compartmentNo];
    CellMLAdapterType &cellmlAdapter = heunScheme.discretizableInTime();

    // loop over CellML states
    for (int stateNo = 0; stateNo < nStates; stateNo++)
    {
      //cellmlAdapter.data().states()->setRepresentationGlobal();

      // get the values from the input vector `data` to localValues
      for (dof_no_t dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts; dofNoLocal++)
      {
        localValues[dofNoLocal] = data[dofNoLocal * nCompartments * nStates + compartmentNo * nStates + stateNo];
      }

      // store all the values for the current state of the current compartment
      cellmlAdapter.data().states()->setValuesWithoutGhosts(stateNo, localValues);
    }
  }
}

template<typename StrangSplittingMultidomain>
void MultidomainWrapper<StrangSplittingMultidomain>::
reset()
{
  strangSplittingMultidomain_.reset();
}

template<typename StrangSplittingMultidomain>
typename MultidomainWrapper<StrangSplittingMultidomain>::Data &MultidomainWrapper<StrangSplittingMultidomain>::
data()
{
  // get a reference to the data object
  return strangSplittingMultidomain_.data();
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class
template<typename StrangSplittingMultidomain>
std::shared_ptr<typename MultidomainWrapper<StrangSplittingMultidomain>::OutputConnectorDataType> MultidomainWrapper<StrangSplittingMultidomain>::
getOutputConnectorData()
{
  return strangSplittingMultidomain_.getOutputConnectorData();
}


//! set a new time step width, gets transferred to numberTimeSteps_
template<typename StrangSplittingMultidomain>
void MultidomainWrapper<StrangSplittingMultidomain>::
setTimeStepWidth(double timeStepWidth)
{
  return strangSplittingMultidomain_.setTimeStepWidth(timeStepWidth);
}

//! set a new number of time steps
template<typename StrangSplittingMultidomain>
void MultidomainWrapper<StrangSplittingMultidomain>::
setNumberTimeSteps(int numberTimeSteps)
{
  return strangSplittingMultidomain_.setNumberTimeSteps(numberTimeSteps);
}

//! set a new time interval that will be simulated by next call to advanceTimeSpan. This also potentially changes the time step width (it preserves the number of timesteps in the new time span)
template<typename StrangSplittingMultidomain>
void MultidomainWrapper<StrangSplittingMultidomain>::
setTimeSpan(double startTime, double endTime)
{
  return strangSplittingMultidomain_.setTimeSpan(startTime, endTime);
}

//! interval for output of time step number and time
template<typename StrangSplittingMultidomain>
int MultidomainWrapper<StrangSplittingMultidomain>::
timeStepOutputInterval()
{
  return strangSplittingMultidomain_.timeStepOutputInterval();
}
//! start time of time interval to be simulated
template<typename StrangSplittingMultidomain>
double MultidomainWrapper<StrangSplittingMultidomain>::
startTime()
{
  return strangSplittingMultidomain_.startTime();
}

//! end time of simulation
template<typename StrangSplittingMultidomain>
double MultidomainWrapper<StrangSplittingMultidomain>::
endTime()
{
  return strangSplittingMultidomain_.endTime();
}

//! number of time steps in simulation time
template<typename StrangSplittingMultidomain>
int MultidomainWrapper<StrangSplittingMultidomain>::
numberTimeSteps()
{
  return strangSplittingMultidomain_.numberTimeSteps();
}

//! time step for simulation
template<typename StrangSplittingMultidomain>
double MultidomainWrapper<StrangSplittingMultidomain>::
timeStepWidth()
{
  return strangSplittingMultidomain_.timeStepWidth();
}

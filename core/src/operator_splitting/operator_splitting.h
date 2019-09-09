#pragma once

#include "time_stepping_scheme/time_stepping_scheme.h"
#include "output_writer/manager.h"
#include "interfaces/runnable.h"
#include "data_management/time_stepping/time_stepping.h"
#include "partition/rank_subset.h"
#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
class OperatorSplitting :
  public ::TimeSteppingScheme::TimeSteppingScheme,    // contains also Multipliable
  public Runnable
  //public Printer<typename TimeStepping2::OutputConnectorDataType>
{
public:
  typedef typename TimeStepping1::FunctionSpace FunctionSpace;
  typedef typename TimeStepping1::Data Data;
  typedef typename TimeStepping1::OutputConnectorDataType OutputConnectorDataType;  // needed when this class is itself part of an operator splitting
 
  //! constructor
  OperatorSplitting(DihuContext context, std::string schemeName);

  //! destructor
  virtual ~OperatorSplitting() {}

  //! run the simulation
  void run();

  //! get the data to be reused in further computations
  OutputConnectorDataType getOutputConnectorData();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);
  
  //! initialize data
  void initialize();

  //! reset state such that new initialization becomes necessary
  virtual void reset();

  //! return the data object
  Data &data();

  //! get a reference to the first timestepping object
  TimeStepping1 &timeStepping1();

  //! get a reference to the second timestepping object
  TimeStepping2 &timeStepping2();

  //! output the given data for debugging
  std::string getString(OutputConnectorDataType &data);

protected:

  TimeStepping1 timeStepping1_;    ///< the object to be discretized
  TimeStepping2 timeStepping2_;    ///< the object to be discretized

  int timeStepOutputInterval_;    ///< time step number and time is output every timeStepOutputInterval_ time steps
  std::string schemeName_;        ///< the key as in the contig, i.e. "Strang" or "Godunov" or "Coupling", only for debugging outputs
  std::string logKeyTimeStepping1AdvanceTimeSpan_;  ///< key for logging of the duration of the advanceTimeSpan() call of timeStepping1
  std::string logKeyTimeStepping2AdvanceTimeSpan_;  ///< key for logging of the duration of the advanceTimeSpan() call of timeStepping2
  std::string logKeyTransfer12_;  ///< key for logging of the duration of data transfer from timestepping 1 to 2
  std::string logKeyTransfer21_;  ///< key for logging of the duration of data transfer from timestepping 2 to 1

  std::string transferSlotName_;  ///< some solver objects have multiple output slots, e.g. cellMLAdapter has intermediates and states as possible output values to use for further computation. transferSlotName select which one to use in the transfer of this operator splitting

  bool initialized_;               ///< if initialize() was already called
};

/*
template<typename OutputConnectorDataType>
class Printer
{
  void print(OutputConnectorDataType &data);
};

*/

}  // namespace

#include "operator_splitting/operator_splitting.tpp"

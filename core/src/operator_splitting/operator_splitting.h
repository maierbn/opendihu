#pragma once

#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "output_writer/manager.h"
#include "interfaces/runnable.h"
#include "data_management/time_stepping/time_stepping.h"
#include "partition/rank_subset.h"
#include "output_connector_data_transfer/output_connector_data_transfer.h"
#include "data_management/operator_splitting.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
class OperatorSplitting :
  public ::TimeSteppingScheme::TimeSteppingScheme,    // contains also Multipliable
  public Runnable
{
public:
  typedef typename TimeStepping1::FunctionSpace FunctionSpace;
  typedef Data::OperatorSplitting<TimeStepping1,TimeStepping2> Data;
  typedef typename Data::OutputConnectorDataType OutputConnectorDataType;  // needed when this class is itself part of an operator splitting
  typedef TimeStepping1 TimeStepping1Type;
  typedef TimeStepping2 TimeStepping2Type;

  //! constructor
  OperatorSplitting(DihuContext context, std::string schemeName);

  //! destructor
  virtual ~OperatorSplitting() {}

  //! run the simulation
  void run();

  //! get the data to be reused in further computations
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

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
  std::string getString(std::shared_ptr<OutputConnectorDataType> data);

protected:

  TimeStepping1 timeStepping1_;     //< the object to be discretized
  TimeStepping2 timeStepping2_;     //< the object to be discretized

  Data data_;                       //< data object that stores the outputConnectorData_ object which is a tuple of both outputConnectorData objects of the timestepping schemes

  int timeStepOutputInterval_;      //< time step number and time is output every timeStepOutputInterval_ time steps
  std::string schemeName_;          //< the key as in the contig, i.e. "Strang" or "Godunov" or "Coupling", only for debugging outputs
  std::string description_;         //< a description that will be printed in debugging output and in the solver structure visualization
  std::string logKeyTimeStepping1AdvanceTimeSpan_;  //< key for logging of the duration of the advanceTimeSpan() call of timeStepping1
  std::string logKeyTimeStepping2AdvanceTimeSpan_;  //< key for logging of the duration of the advanceTimeSpan() call of timeStepping2
  std::string logKeyTransfer12_;    //< key for logging of the duration of data transfer from timestepping 1 to 2
  std::string logKeyTransfer21_;    //< key for logging of the duration of data transfer from timestepping 2 to 1

  std::shared_ptr<OutputConnection> outputConnection_; //< information regarding the mapping between the data slots of the two terms

  bool initialized_;                //< if initialize() was already called
};

}  // namespace

#include "operator_splitting/operator_splitting.tpp"

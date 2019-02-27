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
  //public Printer<typename TimeStepping2::TransferableSolutionDataType>
{
public:
  typedef typename TimeStepping1::FunctionSpace FunctionSpace;
  typedef typename TimeStepping1::Data Data;
  typedef typename TimeStepping1::TransferableSolutionDataType TransferableSolutionDataType;  // needed when this class is itself part of an operator splitting
 
  //! constructor
  OperatorSplitting(DihuContext context, std::string schemeName);

  //! destructor
  virtual ~OperatorSplitting() {}

  //! run the simulation
  void run();

  //! get the data to be reused in further computations
  TransferableSolutionDataType getSolutionForTransfer();

  //! return whether the object has a specified mesh type or if it is independent of any mesh type
  bool knowsMeshType();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);
  
  //! initialize data
  void initialize();

  //! reset state such that new initialization becomes necessary
  virtual void reset();

  //! return the data object
  Data &data();

  //! output the given data for debugging
  std::string getString(TransferableSolutionDataType &data);

protected:

  TimeStepping1 timeStepping1_;    ///< the object to be discretized
  TimeStepping2 timeStepping2_;    ///< the object to be discretized

  int timeStepOutputInterval_;    ///< time step number and time is output every timeStepOutputInterval_ time steps
  std::string schemeName_;        ///< the key as in the contig, i.e. "Strang" or "Godunov" or "Coupling", only for debugging outputs

  bool initialized_;               ///< if initialize() was already called
};

/*
template<typename TransferableSolutionDataType>
class Printer
{
  void print(TransferableSolutionDataType &data);
};

*/

}  // namespace

#include "operator_splitting/operator_splitting.tpp"

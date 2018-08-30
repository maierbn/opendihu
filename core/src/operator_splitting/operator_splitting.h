#pragma once

#include "time_stepping_scheme/time_stepping_scheme.h"
#include "output_writer/manager.h"
#include "control/runnable.h"
#include "data_management/time_stepping.h"
#include "partition/rank_subset.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
class OperatorSplitting :
  public ::TimeSteppingScheme::TimeSteppingScheme, public Runnable
{
public:
  typedef typename TimeStepping1::FunctionSpace FunctionSpace;
  typedef typename TimeStepping1::Data Data;
 
  //! constructor
  OperatorSplitting(DihuContext context, std::string schemeName);

  //! destructor
  virtual ~OperatorSplitting() {}

  //! run the simulation
  void run();

  //! return a solution vector
  Vec &solution();

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

protected:

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer


  TimeStepping1 timeStepping1_;    ///< the object to be discretized
  TimeStepping2 timeStepping2_;    ///< the object to be discretized

  bool outputData1_;               ///< if data output via writer is enabled for timeStepping1
  bool outputData2_;               ///< if data output via writer is enabled for timeStepping2

  bool initialized_;               ///< if initialize() was already called
};

}  // namespace

#include "operator_splitting/operator_splitting.tpp"
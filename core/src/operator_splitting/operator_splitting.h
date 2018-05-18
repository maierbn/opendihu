#pragma once

#include "time_stepping_scheme/time_stepping_scheme.h"
#include "output_writer/manager.h"
#include "control/runnable.h"
#include "data_management/time_stepping.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
class OperatorSplitting :
  public ::TimeSteppingScheme::TimeSteppingScheme, public Runnable
{
public:
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


protected:

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer

  void initialize();

  TimeStepping1 timeStepping1_;    ///< the object to be discretized
  TimeStepping2 timeStepping2_;    ///< the object to be discretized

  bool outputData1_;               ///< if data output via writer is enabled for timeStepping1
  bool outputData2_;               ///< if data output via writer is enabled for timeStepping2

};

}  // namespace

#include "operator_splitting/operator_splitting.tpp"
#pragma once

#include <Python.h>  // has to be the first included header

#include "control/dihu_context.h"
#include "control/python_config/python_config.h"

namespace Control
{

/** This class stores the times when different fibers are stimulated and output them to a log file at the end.
  */
class StimulationLogging
{
public:

  //! constructor
  StimulationLogging(PythonConfig specificSettings);

  //! register the event of starting a stimulation for given motor unit, called by FastMonodomainSolver
  static void logStimulationBegin(double currentTime, int motorUnitNo, int fiberNo);

  //! this will be called at the end of the simulation run
  static void writeLogFile();

  struct StimulationLogEntry
  {
    double time;      //< time when stimulation starts
    int motorUnitNo;  //< motor unit number that is stimulated, set to -1 if not set
    int fiberNo;      //< fiber number , set to -1 if not set
  };

private:

  static std::string filename_;    //< filename of the log file

  static std::vector<StimulationLogEntry> logEntries_;      //< all entries that will be written from the current MPI rank to the log file
};

//! output operator for StimulationLogEntry
std::ostream &operator<<(std::ostream &stream, const StimulationLogging::StimulationLogEntry rhs);

}  // namespace

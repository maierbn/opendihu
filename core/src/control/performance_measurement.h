#pragma once

#include <Python.h>  // has to be the first included header
#include <map>

#include "dihu_context.h"
#include "interfaces/runnable.h"

namespace Control
{

/** A class used for timing and error performance measurements. Timing is done using MPI_Wtime.
 */
class PerformanceMeasurement
{
public:

  //! start timing measurement for a given keyword
  static void start(std::string name);

  //! stop timing measurement for a given keyword, the counter of number of time spans is increased by numberAccumulated
  static void stop(std::string name, int numberAccumulated=1);

  //! compute the mean magnitude of the given error vector or matrix and store it under name
  template<typename T>
  static void measureError(std::string name, T differenceVector);

  //! write collected information to a log file
  static void writeLogFile(std::string logFileName = "logs/log");

  //! save a parameter that will be added in the log file
  template<typename T>
  static void setParameter(std::string key, T parameter);

private:

  //! parse some system information
  static void parseStatusInformation();

  struct Measurement
  {
    //! constructor
    Measurement();

    double start;   ///< last start time point
    double totalDuration;   ///< sum of previous measurements
    int nTimeSpans;     ///< the number of measurements that lead to the total time in totalDuration

    double totalError;  ///< sum of all errors
    int nErrors;        ///< number of summands of totalError
  };

  static std::map<std::string, Measurement> measurements_;   ///< the currently stored measurements
  static std::map<std::string,std::string> parameters_;   ///< arbitrary parameters that will be stored in the log


};

template<>
void PerformanceMeasurement::measureError<double>(std::string name, double differenceVector);

};  // namespace
#include "control/performance_measurement.tpp"

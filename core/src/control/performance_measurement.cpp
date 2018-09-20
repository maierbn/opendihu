#include "control/performance_measurement.h"

#include "easylogging++.h"
#include <memory>
#include <fstream>

namespace Control
{

std::map<std::string, PerformanceMeasurement::Measurement> PerformanceMeasurement::measurements_;

PerformanceMeasurement::Measurement::Measurement() :
  start(0.0), totalDuration(0.0), nTimeSpans(0), totalError(0.0), nErrors(0)
{
}

void PerformanceMeasurement::start(std::string name)
{
  std::map<std::string, Measurement>::iterator iter = measurements_.find(name);

  // if there is no entry of name yet, create new
  if (iter == measurements_.end())
  {
    auto insertedIter = measurements_.insert(std::pair<std::string, Measurement>(name, Measurement()));
    iter = insertedIter.first;
  }

  // measure current time
  iter->second.start = MPI_Wtime();
}

void PerformanceMeasurement::stop(std::string name, int numberAccumulated)
{
  double stopTime = MPI_Wtime();

  if(measurements_.find(name) == measurements_.end())
  {
    LOG(ERROR) << "PerformanceMeasurement stop with name \"" << name << "\", a corresponding start is not present.";
  }
  else
  {
    Measurement &measurement = measurements_[name];
    double duration = stopTime - measurement.start;
    measurement.totalDuration += duration;
    measurement.nTimeSpans += numberAccumulated;
  }
}

void PerformanceMeasurement::log(std::string logFileName)
{
  // determine file name
  std::stringstream filename;
  filename << logFileName << "." << std::setw(7) << std::setfill('0') << DihuContext::ownRankNo() << ".csv";
  logFileName = filename.str();

  // open log file
  std::ofstream file;
  file.open(logFileName, std::ios::out | std::ios::binary | std::ios::trunc);

  if (!file.is_open())
  {
    LOG(ERROR) << "Could not open log file \"" << logFileName << "\" for writing";
    return;
  }

  // write datasets to file
  file << "#name;duration;error;number time spans;number error measurements" << std::endl;

  for (auto &measurement : measurements_)
  {
    file << measurement.first << ";"
     << measurement.second.totalDuration / measurement.second.nTimeSpans << ";"
     << measurement.second.totalError / measurement.second.nErrors << ";"
     << measurement.second.nTimeSpans << ";"
     << measurement.second.nErrors << std::endl;
  }
  file.close();
}

template<>
void PerformanceMeasurement::measureError<double>(std::string name, double differenceVector)
{
  std::map<std::string, Measurement>::iterator iter = measurements_.find(name);

  // if there is no entry of name yet, create new
  if (iter == measurements_.end())
  {
    auto insertedIter = measurements_.insert(std::pair<std::string, Measurement>(name, Measurement()));
    iter = insertedIter.first;
  }

  iter->second.totalError += fabs(differenceVector);
  iter->second.nErrors++;
}
};  // namespace

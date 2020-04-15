#include "control/diagnostic_tool/performance_measurement.h"

#include "easylogging++.h"
#include <memory>
#include <fstream>
#include <unistd.h>
#include <limits.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/param.h>
#include <iomanip>
//#include <stdlib.h>  //was only for function getenv()

#include "output_writer/generic.h"

namespace Control
{

std::map<std::string, PerformanceMeasurement::Measurement> PerformanceMeasurement::measurements_;
std::map<std::string,std::string> PerformanceMeasurement::parameters_;
std::map<std::string, int> PerformanceMeasurement::sums_;
std::shared_ptr<std::thread> PerformanceMeasurement::perfThread_;
int PerformanceMeasurement::perfThreadHandle_;

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

  if (measurements_.find(name) == measurements_.end())
  {
    LOG(ERROR) << "PerformanceMeasurement stop with name \"" << name << "\", a corresponding start is not present.";
  }
  else
  {
    Measurement &measurement = measurements_[name];
    double duration = stopTime - measurement.start;
    measurement.totalDuration += duration;
    measurement.nTimeSpans += numberAccumulated;

    VLOG(2) << "PerformanceMeasurement::stop(" << name << "), time span [" << measurement.start << "," << stopTime << "], duration=" << duration
      << ", now total: " << measurement.totalDuration << ", nTimeSpans: " << measurement.nTimeSpans;
  }
}

void PerformanceMeasurement::startFlops()
{
  // get process id
  int pid = getpid();

  perfThread_ = std::make_shared<std::thread>([pid](){

    std::stringstream filename;
    filename << "perf." << DihuContext::ownRankNoCommWorld() << ".txt";

    std::stringstream command;
    command << "perf stat -e r5301c7 -p " << pid << " -o " << filename.str();

    LOG(INFO) << "starting perf with command " << command.str();

    int returnValue = system(command.str().c_str());
    LOG(DEBUG) << returnValue;
  });
  perfThreadHandle_ = perfThread_->native_handle();

  perfThread_->detach();
}

void PerformanceMeasurement::endFlops()
{
  pthread_cancel(perfThreadHandle_);

  std::stringstream filename;
  filename << "perf." << DihuContext::ownRankNoCommWorld() << ".txt";

  // parse file
  std::ifstream file(filename.str());
  if (file.is_open())
  {
     std::string content((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));
     LOG(DEBUG) << content;
  }
}

std::string PerformanceMeasurement::getParameter(std::string key)
{
  if (parameters_.find(key) == parameters_.end())
    return "";

  return parameters_[key];
}

double PerformanceMeasurement::getDuration(std::string measurementName)
{
  if (measurements_.find(measurementName) == measurements_.end())
    return 0.0;

  return measurements_[measurementName].totalDuration;
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

void PerformanceMeasurement::countNumber(std::string name, int number)
{
  std::map<std::string, int>::iterator iter = sums_.find(name);
  
  // if there is no entry of name yet, create new
  if (iter == sums_.end())
  {
    auto insertedIter = sums_.insert(std::pair<std::string, int>(name, 0));
    iter = insertedIter.first;
  }
  
  iter->second += number;
}

void PerformanceMeasurement::getMemoryConsumption(int &pageSize, long long &virtualMemorySize, long long &residentSetSize, long long &dataSize, double &totalUserTime)
{
  // adapted from https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c

  std::string pid, comm, state, ppid, pgrp, session, ttyNr, tpgid, flags, minflt, cminflt, majflt, cmajflt, stime, cutime, cstime,
    priority, nice, o, itrealvalue, starttime, size, resident, shared, text, none;

  std::ifstream stream("/proc/self/stat", std::ios_base::in);

  unsigned long vsize;
  long long rss;
  long long utime;
  long long data;

  // parse entries in /proc/self/stat
  stream >> pid >> comm >> state >> ppid >> pgrp >> session >> ttyNr
              >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
              >> utime >> stime >> cutime >> cstime >> priority >> nice
              >> o >> itrealvalue >> starttime >> vsize >> rss;

  stream.close();
  stream.open("/proc/self/statm", std::ios_base::in);

  // parse entries in /proc/self/statm
  stream >> size >> resident >> shared >> text >> none >> data;
  stream.close();

  // compute final values
  pageSize = getpagesize(); // page size in KB
  virtualMemorySize = vsize;
  residentSetSize = rss * pageSize;
  dataSize = data * pageSize;

  totalUserTime = (double)utime / sysconf(_SC_CLK_TCK);

  // http://man7.org/linux/man-pages/man5/proc.5.html
  /** Resident Set Size: number of pages the process has
                      in real memory.  This is just the pages which count
                      toward text, data, or stack space.  This does not
                      include pages which have not been demand-loaded in,
                      or which are swapped out.
  */

}

void PerformanceMeasurement::parseStatusInformation()
{
  // get current memory consumption
  int pageSize;
  long long virtualMemorySize;
  long long residentSetSize;
  long long dataSize;
  double totalUsertime;

  getMemoryConsumption(pageSize, virtualMemorySize, residentSetSize, dataSize, totalUsertime);

  // store parameters
  PerformanceMeasurement::setParameter("memoryPage", pageSize);
  PerformanceMeasurement::setParameter("memoryVirtual", virtualMemorySize);
  PerformanceMeasurement::setParameter("memoryResidentSet", residentSetSize);
  PerformanceMeasurement::setParameter("memoryData", dataSize);
  PerformanceMeasurement::setParameter("totalUsertime", totalUsertime);

  // format total user time in readable format for output
  int minTotal = int(totalUsertime/60);
  int hTotal = int(totalUsertime/3600);
  int d = int(totalUsertime/(3600*24));
  int h = hTotal - d*24;
  int min = minTotal - h*60 - d*24*60;
  double sDouble = totalUsertime - minTotal*60;
  int s = int(sDouble);
  sDouble -= s;

  std::stringstream message;
  if (d != 0)
    message << d << "d ";
  if (h != 0)
    message << h << ":";
  if (min != 0)
  {
    std::stringstream secondsFraction;
    secondsFraction << std::setprecision(4) << sDouble;
    if (h == 0)
      message << min;
    else
      message << std::setw(2) << std::setfill('0') << min;

    message << ":" << std::setw(2) << std::setfill('0') << s << secondsFraction.str().substr(1);
    if (h == 0)
      message << " min";
  }
  else
  {
    message << s << "s";
  }
  LOG(INFO) << "Total user time: " << message.str();
}

} // namespace

#include "control/diagnostic_tool/memory_leak_finder.h"

#include "control/diagnostic_tool/performance_measurement.h"

namespace Control
{

long long int MemoryLeakFinder::currentMemoryConsumption_ = 0;   //< the current number of bytes allocated in residual set memory

long long int MemoryLeakFinder::nBytesIncreaseSinceLastCheck()
{
  // get current memory consumption
  int pageSize;
  long long virtualMemorySize;
  long long residentSetSize;
  long long dataSize;
  double totalUserTime;
  PerformanceMeasurement::getMemoryConsumption(pageSize, virtualMemorySize, residentSetSize, dataSize, totalUserTime);

  long long increment = residentSetSize - currentMemoryConsumption_;
  currentMemoryConsumption_ = residentSetSize;

  return increment;
}

void MemoryLeakFinder::warnIfMemoryConsumptionIncreases(std::string message)
{
  bool firstCall = currentMemoryConsumption_ == 0;
  long long increase = nBytesIncreaseSinceLastCheck();

  if (increase > 1024*1024 && !firstCall)
  {
    LOG(WARNING) << message << ": Memory consumption increased by " << increase << " B to " << currentMemoryConsumption_ << " B.";
  }
}

}  // namespace

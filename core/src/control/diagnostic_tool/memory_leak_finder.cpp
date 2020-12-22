#include "control/diagnostic_tool/memory_leak_finder.h"

#include "control/diagnostic_tool/performance_measurement.h"
#include <iomanip>

namespace Control
{

long long int MemoryLeakFinder::currentMemoryConsumption_ = 0;   //< the current number of kilobytes allocated in residual set memory

long long int MemoryLeakFinder::currentMemoryConsumptionKiloBytes()
{
  // get current memory consumption
  int pageSize;
  long long virtualMemorySize;
  long long residentSetSize;
  long long dataSize;
  double totalUserTime;
  PerformanceMeasurement::getMemoryConsumption(pageSize, virtualMemorySize, residentSetSize, dataSize, totalUserTime);

  return residentSetSize / 1024;
}

long long int MemoryLeakFinder::nKiloBytesIncreaseSinceLastCheck()
{
  // get current memory consumption
  long long residentSetSize = currentMemoryConsumptionKiloBytes();
  long long increment = residentSetSize - currentMemoryConsumption_;
  currentMemoryConsumption_ = residentSetSize;

  return increment;
}

void MemoryLeakFinder::warnIfMemoryConsumptionIncreases(std::string message)
{
  bool firstCall = currentMemoryConsumption_ == 0;
  long long increase = nKiloBytesIncreaseSinceLastCheck();

  if (increase >= 1024 && !firstCall)
  {
    LOG(WARNING) << message << ": Memory consumption increased by " << increase/(1024) << " MB to "
      << currentMemoryConsumption_/(1024) << " MB.";
  }
}

}  // namespace

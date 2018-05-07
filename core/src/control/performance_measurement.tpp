#include "control/performance_measurement.h"

namespace Control
{
  
template<typename T>
void PerformanceMeasurement::measureError(std::string name, T differenceVector)
{
  for (int i = 0; i < differenceVector.size(); i++)
  {
    measureError(name, differenceVector[i]);
  }
}

};  // namespace
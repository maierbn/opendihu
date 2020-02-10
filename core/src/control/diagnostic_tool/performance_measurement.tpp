#include "control/diagnostic_tool/performance_measurement.h"

#include <sstream>

#include "utility/string_utility.h"
#include "utility/vector_operators.h"

namespace Control
{

template<typename T>
void PerformanceMeasurement::
measureError(std::string name, T differenceVector)
{
  for (int i = 0; i < differenceVector.size(); i++)
  {
    measureError(name, differenceVector[i]);
  }
}

template<typename T>
void PerformanceMeasurement::
setParameter(std::string key, T parameter)
{
  std::stringstream str;
  str << parameter;
  parameters_[key] = str.str();
}

} // namespace

#include "time_stepping_scheme/repeated_call_static.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename Solver>
RepeatedCallStatic<Solver>::RepeatedCallStatic(DihuContext context) :
  TimeSteppingScheme(context["RepeatedCallStatic"]), solver_(context_)
{
  this->specificSettings_ = this->context_.getPythonConfig();
}

template<typename Solver>
void RepeatedCallStatic<Solver>::
initialize()
{
  if (initialized_)
    return;
 
  TimeSteppingScheme::initialize();
  LOG(TRACE) << "RepeatedCallStatic::initialize";

  // initialize underlying Solver object, also with time step width
  solver_.initialize();
}

template<typename Solver>
void RepeatedCallStatic<Solver>::advanceTimeSpan()
{
  // run the solver which does not have advanceTimeSpan
  solver_.run();
}

template<typename Solver>
void RepeatedCallStatic<Solver>::
run()
{
  initialize();
  advanceTimeSpan();
}

template<typename Solver>
typename RepeatedCallStatic<Solver>::Data &RepeatedCallStatic<Solver>::
data()
{
  return solver_.data();
}

template<typename Solver>
std::shared_ptr<typename RepeatedCallStatic<Solver>::OutputConnectorDataType> RepeatedCallStatic<Solver>::
getOutputConnectorData()
{
  return solver_.getOutputConnectorData();
}

} // namespace

#include "time_stepping_scheme.h"

#include "control/python_utility.h"

namespace TimeSteppingScheme
{
  
TimeSteppingScheme::TimeSteppingScheme(const DihuContext &context) : 
  context_(context)
{
  specificSettings_ = NULL;   // needs to be set by deriving class
}

void TimeSteppingScheme::setTimeStepWidth(double timeStepWidth)
{
  double epsilon = 1e-15;
  numberTimeSteps_ = int(std::ceil((endTime_ - startTime_) / timeStepWidth - epsilon));
}

void TimeSteppingScheme::setNumberTimeSteps(int numberTimeSteps)
{
  numberTimeSteps_ = numberTimeSteps;
}

void TimeSteppingScheme::setTimeSpan(double startTime, double endTime)
{
  startTime_ = startTime;
  endTime_ = endTime;
}

void TimeSteppingScheme::initialize()
{
  // initialize time stepping values
  startTime_ = 0.0;
  endTime_ = PythonUtility::getOptionDouble(specificSettings_, "endTime", 1.0, PythonUtility::Positive);
  
  if (PythonUtility::containsKey(specificSettings_, "timeStepWidth"))
  {
    setTimeStepWidth(PythonUtility::getOptionDouble(specificSettings_, "timeStepWidth", 0.001, PythonUtility::Positive));
    
    if (PythonUtility::containsKey(specificSettings_, "numberTimeSteps"))
    {
      numberTimeSteps_ = PythonUtility::getOptionInt(specificSettings_, "numberTimeSteps", 10, PythonUtility::Positive);
      LOG(WARNING) << "Time step width will be overridden by number of time steps (" << numberTimeSteps_ << ")";
    }
  }
  else
  {
    numberTimeSteps_ = PythonUtility::getOptionInt(specificSettings_, "numberTimeSteps", 10, PythonUtility::Positive);
  }
  
  LOG(INFO) << "Time span: [" << startTime_ << "," << endTime_ << "], Number of time steps: " << numberTimeSteps_
    << ", time step width: " << (endTime_ - startTime_) / numberTimeSteps_;
}

};  // namespace
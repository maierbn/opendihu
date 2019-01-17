#include "time_stepping_scheme.h"

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

TimeSteppingScheme::TimeSteppingScheme(DihuContext context) :
  Splittable(), context_(context), specificSettings_(NULL), initialized_(false)
{
  // specificSettings_ needs to be set by deriving class, in time_stepping_scheme_ode.tpp
  isTimeStepWidthSignificant_ = false;
}

void TimeSteppingScheme::setTimeStepWidth(double timeStepWidth)
{
  double epsilon = 1e-15;
  numberTimeSteps_ = int(std::ceil((endTime_ - startTime_) / timeStepWidth - epsilon));
  timeStepWidth_ = timeStepWidth;
}

void TimeSteppingScheme::setNumberTimeSteps(int numberTimeSteps)
{
  numberTimeSteps_ = numberTimeSteps;
  timeStepWidth_ = (endTime_ - startTime_) / numberTimeSteps;
  LOG(DEBUG) << "timeStepWidth_ in setNumberTimeSteps: " << timeStepWidth_;
}

void TimeSteppingScheme::setTimeSpan(double startTime, double endTime)
{
  startTime_ = startTime;
  endTime_ = endTime;

  if (isTimeStepWidthSignificant_)
  {
    setTimeStepWidth(timeStepWidth_);
    LOG(DEBUG) << "set number of time steps to " <<numberTimeSteps_<< " from timeStepWidth " << timeStepWidth_;
  }
}

void TimeSteppingScheme::reset()
{
  initialized_ = false;
}

void TimeSteppingScheme::initialize()
{
  if (initialized_)
    return;
  
  LOG(DEBUG) << "TimeSteppingScheme::initialize()";
  
  // initialize time stepping values
  startTime_ = 0.0;
  endTime_ = 1.0;
  if (specificSettings_.hasKey("endTime"))
    endTime_ = specificSettings_.getOptionDouble("endTime", 1.0, PythonUtility::Positive);

  LOG(DEBUG) << "  TimeSteppingScheme::initialize read endTime=" << endTime_;

  if (specificSettings_.hasKey("timeStepWidth"))
  {
    timeStepWidth_ = specificSettings_.getOptionDouble("timeStepWidth", 0.001, PythonUtility::Positive);
    setTimeStepWidth(timeStepWidth_);

    LOG(DEBUG) << "  TimeSteppingScheme::initialize, timeStepWidth="
      << specificSettings_.getOptionDouble("timeStepWidth", 0.001, PythonUtility::Positive)
      << ", compute numberTimeSteps=" <<numberTimeSteps_;

    if (specificSettings_.hasKey("numberTimeSteps"))
    {
      numberTimeSteps_ = specificSettings_.getOptionInt("numberTimeSteps", 10, PythonUtility::Positive);
      isTimeStepWidthSignificant_ = false;
      LOG(WARNING) << "Time step width will be overridden by number of time steps (" << numberTimeSteps_ << ")";

      setNumberTimeSteps(numberTimeSteps_);
    }
    else
    {
      isTimeStepWidthSignificant_ = true;
    }
  }
  else
  {
    int numberTimeSteps = specificSettings_.getOptionInt("numberTimeSteps", 10, PythonUtility::Positive);
    LOG(DEBUG) << "  TimeSteppingScheme::initialize, timeStepWidth not specified, read numberTimeSteps: " << numberTimeSteps;
    setNumberTimeSteps(numberTimeSteps);
  }

  LOG(INFO) << "Time span: [" << startTime_ << "," << endTime_ << "], Number of time steps: " << numberTimeSteps_
    << ", time step width: " << timeStepWidth_;

  // log timeStepWidth as the key that is given by "logTimeStepWidthAsKey"
  if (specificSettings_.hasKey("logTimeStepWidthAsKey"))
  {
    std::string timeStepWidthKey = specificSettings_.getOptionString("logTimeStepWidthAsKey", "timeStepWidth");
    Control::PerformanceMeasurement::setParameter(timeStepWidthKey, timeStepWidth_);
  }

  if (specificSettings_.hasKey("logTimeStepWidthAsKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }

  timeStepOutputInterval_ = specificSettings_.getOptionInt("timeStepOutputInterval", 100, PythonUtility::Positive);

  initialized_ = true;
}

int TimeSteppingScheme::timeStepOutputInterval()
{
  return this->timeStepOutputInterval_;
}

double TimeSteppingScheme::startTime()
{
  return startTime_;
}

double TimeSteppingScheme::endTime()
{
  return endTime_;
}

double TimeSteppingScheme::numberTimeSteps()
{
  return numberTimeSteps_;
}
  
double TimeSteppingScheme::timeStepWidth()
{
  return timeStepWidth_;
}

PythonConfig TimeSteppingScheme::specificSettings()
{
  return specificSettings_;
}

} // namespace

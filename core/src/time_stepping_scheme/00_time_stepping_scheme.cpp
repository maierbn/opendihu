#include "time_stepping_scheme/00_time_stepping_scheme.h"

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
  double epsilon = 1e-1;
  // Increase time step width by maximum of epsilon=10%.
  // Increasing time step width is potentially dangerous, because it can make the timestepping scheme unstable.
  // Ideally changing the timestep width should not be necessary at all, if it is given correctly in the config.
  // Examples:
  //  t_end = 1.09, dt = 1.0 -> increase dt to 1.09, one timestep
  //  t_end = 1.11, dt = 1.0 -> decrease dt to 0.555, two timesteps

  long double n = std::ceil((long double)(endTime_ - startTime_) / timeStepWidth - epsilon);
  numberTimeSteps_ = std::max(1,int(n));

  if (n > (double)std::numeric_limits<int>::max())
  {
    LOG(FATAL) << "Number of timesteps " << n << " for timeStepWidth " << timeStepWidth
      << " in time span [" << startTime_ << "," << endTime_ << "] is too high.";
  }

  setNumberTimeSteps(numberTimeSteps_);
}

void TimeSteppingScheme::setNumberTimeSteps(int numberTimeSteps)
{
  numberTimeSteps_ = numberTimeSteps;
  timeStepWidth_ = (endTime_ - startTime_) / numberTimeSteps;
  LOG(DEBUG) << "numberTimeSteps_: " << numberTimeSteps_ << ", endTime:" << endTime_
    << ", startTime: " << startTime_ << " timeStepWidth_ in setNumberTimeSteps: " << timeStepWidth_;
}

void TimeSteppingScheme::setTimeSpan(double startTime, double endTime)
{
  startTime_ = startTime;
  endTime_ = endTime;

  if (timeStepWidth_ > endTime_-startTime_)
  {
    LOG(DEBUG) << "time span [" << startTime << "," << endTime << "], reduce timeStepWidth from " << timeStepWidth_ << " to " << endTime_-startTime_;
    timeStepWidth_ = endTime_-startTime_;
  }

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

  VLOG(1) << "in TimeSteppingScheme::initalize, specificSettings_: " << specificSettings_;
  if (specificSettings_.hasKey("endTime"))
    endTime_ = specificSettings_.getOptionDouble("endTime", 1.0, PythonUtility::Positive);

  LOG(DEBUG) << "  TimeSteppingScheme::initialize read endTime=" << endTime_;

  if (specificSettings_.hasKey("timeStepWidth"))
  {
    timeStepWidth_ = specificSettings_.getOptionDouble("timeStepWidth", 0.001, PythonUtility::Positive);
    setTimeStepWidth(timeStepWidth_);

    LOG(DEBUG) << "  TimeSteppingScheme::initialize, timeStepWidth in settings: "
      << specificSettings_.getOptionDouble("timeStepWidth", 0.001, PythonUtility::Positive)
      << ", timeStepWidth_: " << timeStepWidth_ << ", compute numberTimeSteps=" << numberTimeSteps_;

    if (specificSettings_.hasKey("numberTimeSteps"))
    {
      const int numberTimeSteps = specificSettings_.getOptionInt("numberTimeSteps", 10, PythonUtility::NonNegative);
      if (numberTimeSteps > 0)
      {
        numberTimeSteps_ = numberTimeSteps;
        isTimeStepWidthSignificant_ = false;
        LOG(WARNING) << "Time step width (" << timeStepWidth_ << ") will be overridden by number of time steps (" << numberTimeSteps_ << ")";

        setNumberTimeSteps(numberTimeSteps_);
      }
      else
      {
        LOG(WARNING) << "Time step width (" << timeStepWidth_ << ") will NOT be overridden by number of time steps (" << numberTimeSteps << ")";
        isTimeStepWidthSignificant_ = true;
      }
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

  LOG(DEBUG) << "Time span: [" << startTime_ << "," << endTime_ << "], Number of time steps: " << numberTimeSteps_
    << ", time step width: " << timeStepWidth_;

  // log timeStepWidth as the key that is given by "logTimeStepWidthAsKey"
  if (specificSettings_.hasKey("logTimeStepWidthAsKey"))
  {
    std::string timeStepWidthKey = specificSettings_.getOptionString("logTimeStepWidthAsKey", "timeStepWidth");
    Control::PerformanceMeasurement::setParameter(timeStepWidthKey, timeStepWidth_);
  }

  if (specificSettings_.hasKey("durationLogKey"))
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

int TimeSteppingScheme::numberTimeSteps()
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

OutputWriter::Manager TimeSteppingScheme::outputWriterManager()
{
  return outputWriterManager_;
}

std::string TimeSteppingScheme::durationLogKey()
{
  return durationLogKey_;
}

}  // namespace


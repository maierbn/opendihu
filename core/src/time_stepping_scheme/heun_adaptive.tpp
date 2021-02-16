#include "time_stepping_scheme/heun_adaptive.h"

#include <Python.h>
#include <memory>
#include <iostream>
#include <fstream>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
HeunAdaptive<DiscretizableInTime>::HeunAdaptive(DihuContext context) :
  TimeSteppingExplicit<DiscretizableInTime>(context, "HeunAdaptive")
{
}

template<typename DiscretizableInTime>
void HeunAdaptive<DiscretizableInTime>::initialize()
{
  // create data object for heun
  this->data_ = std::make_shared<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>>(this->context_);

  // initialize already writes the first output file
  TimeSteppingSchemeOde<DiscretizableInTime>::initialize();

  // read tolerance from settings
  tolerance_ = this->specificSettings_.getOptionDouble("tolerance", 0.1, PythonUtility::Positive);

  // initialize delta. Initialized very low, later dynamically set to 10% of timeStepWidth
  delta_ = 1e-6;

  // saved timeStepWidth. Dynamically set in simulation
  savedTimeStepWidth_ = -1;

  // check to determine if minTimeStepWidth is used
  minimum_check_ = false;

  // check to determine if delta is used to advance timestepWidth
  ten_percent_check_ = false;

  // read minimal TimeStepWidth from settings
  minTimeStepWidth_ = this->specificSettings_.getOptionDouble("minTimeStepWidth", 1e-6, PythonUtility::Positive);

  // initiate counter for successive steps
  stepCounterSuccess_ = 0;

  // initiate counter for failed steps
  stepCounterFail_ = 0;

  // get strategy to adapt timeStepWidth. Default value is method with non-uniform steps
  timeStepAdaptOption_ = this->specificSettings_.getOptionString("timeStepAdaptOption", "regular");

  // read lowest multiplier from settings
  lowestMultiplier_ = this->specificSettings_.getOptionInt("lowestMultiplier", 1000, PythonUtility::Positive);

  // initialize passed time
  currentTimeHeun_ = 0;

  // initialize multiplicator for the modified method
  multiplicator_ = 0;

  // get the filename of the log file for timestep widths
  timeStepWidthsLogFilename_ = this->specificSettings_.getOptionString("timeStepWidthsLogFilename", "dt_log.csv");
}

//! destructor to write log file
template<typename DiscretizableInTime>
HeunAdaptive<DiscretizableInTime>::
~HeunAdaptive()
{
  // output log file
  if (!timeStepWidthsLogFilename_.empty() && !timeStepWidths_.empty())
  {
    // write log file
    LOG(INFO) << "Create log file \"" << timeStepWidthsLogFilename_ << "\" with " << timeStepWidths_.size() << " timestep widths.";

    std::ofstream logFile;
    OutputWriter::Generic::openFile(logFile, timeStepWidthsLogFilename_);
    for (int i = 0; i < timeStepWidths_.size(); i++)
    {
      logFile << times_[i] << ";" << timeStepWidths_[i] << ";\n";
    }
    logFile.close();
  }
}

template<typename DiscretizableInTime>
void HeunAdaptive<DiscretizableInTime>::
advanceTimeSpan(bool withOutputWritersEnabled)
{

  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timeSpan of current step
  double timeSpan = this->endTime_ - this->startTime_;

  // pre-allocate buffer for logging out timestep widths
  if (timeStepWidths_.empty())
    timeStepWidths_.reserve(1e5);

  // log info
  LOG(DEBUG) << "HeunAdaptive::advanceTimeSpan, timeSpan=" << timeSpan<< ", initial timeStepWidth=" << this->timeStepWidth_;

  // we need to cast the pointer type to the derived class. Otherwise the additional algebraicIncrement()-method of the class TimeSteppingHeun won't be there:
  std::shared_ptr<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>> dataHeun
    = std::static_pointer_cast<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>>(this->data_);

  // get vectors of all components in struct-of-array order, as needed by CellML (i.e. one long vector with [state0 state0 state0 ... state1 state1...]
  Vec &solution = this->data_->solution()->getValuesContiguous();
  Vec &increment = this->data_->increment()->getValuesContiguous();
  Vec &algebraicIncrement = dataHeun->algebraicIncrement()->getValuesContiguous();

  // duplicate current solutions to create same-structured vectors
  VecDuplicate(solution, &temp_solution_normal);
  VecDuplicate(solution, &temp_solution_tilde);
  VecDuplicate(solution, &temp_solution_tilde_algebraic);
  VecDuplicate(increment, &temp_increment_1);
  VecDuplicate(algebraicIncrement, &temp_increment_2);

  // check for savedTimeStepWidth
  if (savedTimeStepWidth_ > 0)
  {
    this->timeStepWidth_ = std::min(savedTimeStepWidth_, timeSpan);
    if (this->timeStepWidth_ == savedTimeStepWidth_)
    {
      LOG(DEBUG) << "SavedTimeStepWidth loaded. Set timeStepWidth to:  " << this->timeStepWidth_;
    }
    else
    {
      LOG(DEBUG) << "Loaded SavedTimeStepWidth to big for timeSpan. Set timeStepWidth to:  " << this->timeStepWidth_;
    }
  }

  // counting time steps in one loop
  int timeStepNo = 0;

  // calculate current time
  double currentTime = this->startTime_;
  LOG(DEBUG) << "New timeSpan started. Current time: " << currentTime;

  // update current time for returning method
  currentTimeHeun_ = currentTimeHeun_ + timeSpan;

  // regular method chosen in config
  if (timeStepAdaptOption_ == "regular")
  {

    // log info
    LOG(DEBUG) << "Option used: " << timeStepAdaptOption_ ;

    // loop over time steps
    for(double time = 0.0; time < timeSpan;)
    {
      if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
      {
        // log info of current progress in timeSpan
        LOG(INFO) << "HeunAdaptive (regular), time " << time << "/" << timeSpan << ", current timeStepWidth: " << this->timeStepWidth_ << ", time remaining: " << timeSpan - time;
      }

      //loop until next solution is found
      while(true)
      {
        // reset vecNorm
        vecnorm_ = 0.0;

        // copy current solution in temporal vectors
        VecCopy(solution, temp_solution_normal);
        VecCopy(solution, temp_solution_tilde);
        VecCopy(solution, temp_solution_tilde_algebraic);

        LOG(DEBUG) << "starting from solution: " << *this->data_->solution();

        // calculate solution with current timeStepWidth. Use modified calculation already used in Heun
        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        temp_solution_normal, temp_increment_1, timeStepNo, currentTime);

        VecAXPY(temp_solution_normal, this->timeStepWidth_, temp_increment_1);

        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        temp_solution_normal, temp_increment_2, timeStepNo + 1, (currentTime + this->timeStepWidth_));

        VecAXPY(temp_increment_2, -1.0, temp_increment_1);
        VecAXPY(temp_solution_normal, 0.5 * this->timeStepWidth_, temp_increment_2);

        // now calculate reference solution consisting of 2 steps with timeStepWidth/2
        // first step
        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        temp_solution_tilde_algebraic, temp_increment_1, timeStepNo, currentTime);

        VecAXPY(temp_solution_tilde_algebraic, 0.5*this->timeStepWidth_, temp_increment_1);

        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        temp_solution_tilde_algebraic, temp_increment_2, timeStepNo + 1, (currentTime + 0.5*this->timeStepWidth_));

        VecAXPY(temp_increment_2, 1.0, temp_increment_1);
        VecAXPY(temp_solution_tilde, 0.25*this->timeStepWidth_, temp_increment_2);

        // second step
        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        temp_solution_tilde, temp_increment_1, timeStepNo, (currentTime + 0.5*this->timeStepWidth_));

        VecCopy(temp_solution_tilde, temp_solution_tilde_algebraic);

        VecAXPY(temp_solution_tilde_algebraic, 0.5*this->timeStepWidth_, temp_increment_1);

        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        temp_solution_tilde_algebraic, temp_increment_2, timeStepNo + 1, (currentTime + this->timeStepWidth_));

        VecAXPY(temp_increment_2, 1.0, temp_increment_1);

        VecAXPY(temp_solution_tilde, 0.25*this->timeStepWidth_, temp_increment_2);

        // check, if solutions are equal
        PetscBool flag = PETSC_FALSE;
        //VecEqual(temp_solution_tilde,temp_solution_normal, &flag);

        // calculate estimator
        if (!flag)
        {
          VecAXPY(temp_solution_tilde, -1.0, temp_solution_normal);
          VecNorm(temp_solution_tilde, NORM_2, &vecnorm_);

          PetscReal vecnorm_Combined = 0;
          MPI_Allreduce(&vecnorm_, &vecnorm_Combined, 1, MPIU_REAL, MPIU_MAX, this->data_->solution()->functionSpace()->meshPartition()->rankSubset()->mpiCommunicator());

          LOG(DEBUG) << "own vecnorm: " << vecnorm_;
          LOG(DEBUG) << ", sum of vecnorms: " << vecnorm_Combined;

          // prevent that the estimator gets nan, if the divisor is 0
          if (fabs(this->timeStepWidth_) < 1e-12)
          {
            estimator_ = 1e-10;
            LOG(DEBUG) << "estimator would be nan, set to 1e-10.";
          }
          else if (fabs(vecnorm_Combined) > 1e+75 || !std::isfinite(vecnorm_Combined))
          {
            LOG(DEBUG) << "vecnorm_Combined (" << vecnorm_Combined << ") is inf, set estimator to 1.2*tolerance_=" << tolerance_*1.2;
            estimator_ = tolerance_ * 1.2;
          }
          else
          {
            estimator_ = (double)vecnorm_Combined / ((1.0 - pow(0.5, 2))*this->timeStepWidth_);
            LOG(DEBUG) << "compute estimator = " << vecnorm_Combined << "/" << ((1.0 - pow(0.5, 2))*this->timeStepWidth_) << " = " << estimator_;
          }
        }
        else
        {
          // Solutions are equal. Set low estimator to make progress
          estimator_ = 1e-10;
          LOG(DEBUG) << "Vecs are equal!";
        }

        // calculate alpha
        alpha_ = pow((tolerance_/estimator_), (1.0/3.0));

        // lof info about calculated values
        LOG(DEBUG) << "With timeStepWidth: " << this->timeStepWidth_ << ", and estimator: " << estimator_ << ", calculated alpha: " << alpha_;

        // bypass calculation and accept step, when one check is set
        if ((minimum_check_ == true) || (ten_percent_check_ == true))
        {

          // reset 10% check
          ten_percent_check_ = false;

          // reset minimum_check only if alpha greater than 1. Prevents rejecting solution after accepted step
          if (alpha_ >= 1.0)
          {
            minimum_check_ = false;
          }
          else
          {
            alpha_ = 1.0;
          }

          // set estimator very low to guarantee progress
          estimator_ = 1e-10;
        }

        // error is below tolerance, accept
        if (estimator_ <= tolerance_)
        {

          // take already calculated value with current timeStepWidth
          VecCopy(temp_solution_normal, solution);

          // apply the prescribed boundary condition values
          this->applyBoundaryConditions();

          // log info
          VLOG(1) << *this->data_->solution();

          // advance simulation time and timeStepNo for the Logging
          timeStepNo++;
          time = time + this->timeStepWidth_;

          // raise counter for successful step
          stepCounterSuccess_ = stepCounterSuccess_ + 1;

          // marker for the output reader
          //LOG(DEBUG) << "XXX" << currentTime << ";" << this->timeStepWidth_;

          // adjust timeStepWidth based on alpha
          this->timeStepWidth_ = 0.9 * alpha_ * this->timeStepWidth_;

          // check if new timeStepWidth is below minimum
          if (this->timeStepWidth_ < minTimeStepWidth_)
          {
            // set check variable
            minimum_check_ = true;
            // adjust timeStepWidth
            this->timeStepWidth_ = minTimeStepWidth_;

            // log info
            LOG(DEBUG) << "New timeStepWidth below minimum. Make next step with " << minTimeStepWidth_;
          }

          // update loop variable
          currentTime = this->startTime_ + time;

          // log info
          LOG(DEBUG) << "Solution accepted. TimeStepWidth set to: " << this->timeStepWidth_;

          // calculate new delta
          delta_ = 0.1 * this->timeStepWidth_;

          // check, if there is still time in timeSpan
          if ((timeSpan - time) != 0.0)
          {
            // remaining time in timeSpan is too low for current timeStep (regarding delta)
            if (((timeSpan - time) - this->timeStepWidth_) <= delta_)
            {
              // timeStepWidth too big, has to be reduced
              if ((timeSpan - time) < this->timeStepWidth_)
              {
                // save current timeStepWidth
                this->savedTimeStepWidth_ = this->timeStepWidth_;

                // set timeStepWidth to remaining time in TimeSpan
                this->timeStepWidth_ = timeSpan - time;

                // log info
                LOG(DEBUG) << "Remaining time for next timeStep to small. Setting timestepwidth to " << this->timeStepWidth_;
              }
              else
              {
                // save current timeStepWidth
                this->savedTimeStepWidth_ = this->timeStepWidth_;

                // set timeStepWidth to remaining time in TimeSpan
                this->timeStepWidth_ = timeSpan - time;

                // set check to accept solution in new step
                ten_percent_check_ = true;

                // log info
                LOG(DEBUG) << "TimeStep short before end of timeSpan. Increase to " << this->timeStepWidth_;
              }
            }
          }
          else
          {
            // log end of timeSpan
            LOG(DEBUG) << "End of timeSpan reached";
          }

          // break inner loop, since next solution was found
          break;
        }
        else    // error is above tolerance, reject
        {
          // adjust timeStepWidth
          this->timeStepWidth_ = 0.9 * alpha_ * this->timeStepWidth_;

          // raise counter for failed step
          stepCounterFail_ = stepCounterFail_ + 1;

          // check if new timeStepWidth is below minimum
          if (this->timeStepWidth_ < minTimeStepWidth_)
          {
            // set check variable
            minimum_check_ = true;

            // adjust timeStepWidth
            this->timeStepWidth_ = minTimeStepWidth_;
          }
          // log info
          LOG(DEBUG) << "Solution rejected. Retry with smaller timeStepWidth " << this->timeStepWidth_;
        }
      }

      // stop duration measurement
      if (this->durationLogKey_ != "")
        Control::PerformanceMeasurement::stop(this->durationLogKey_);

      // log timestep width
      if (!timeStepWidthsLogFilename_.empty())
      {
        timeStepWidths_.push_back(this->timeStepWidth_);
        times_.push_back(currentTime);
      }

      // write current output values
      if (withOutputWritersEnabled)
        this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);

      // start duration measurement
      if (this->durationLogKey_ != "")
        Control::PerformanceMeasurement::start(this->durationLogKey_);
    }
  }
  else if (timeStepAdaptOption_ == "modified") // modified option was chosen
  {

    // log info
    LOG(DEBUG) << "Option used: " << timeStepAdaptOption_ ;

    // get multiplicator for timeStepwidth
    multiplicator_ = floor(timeSpan/this->timeStepWidth_);

    // reset vecNorm
    vecnorm_ = 0.0;

    // copy current solution in temporal vectors
    VecCopy(solution, temp_solution_normal);
    VecCopy(solution, temp_solution_tilde);
    VecCopy(solution, temp_solution_tilde_algebraic);

    VLOG(1) << "starting from solution: " << *this->data_->solution();

    // calculate solution with current timeStepWidth. Use modified calculation already used in Heun
    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    temp_solution_normal, temp_increment_1, timeStepNo, currentTime);

    VecAXPY(temp_solution_normal, this->timeStepWidth_, temp_increment_1);

    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    temp_solution_normal, temp_increment_2, timeStepNo + 1, (currentTime + this->timeStepWidth_));

    VecAXPY(temp_increment_2, -1.0, temp_increment_1);

    VecAXPY(temp_solution_normal, 0.5 * this->timeStepWidth_, temp_increment_2);

    // now calculate reference solution consisting of 2 steps with timeStepWidth/2
    // first step
    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    temp_solution_tilde_algebraic, temp_increment_1, timeStepNo, currentTime);

    VecAXPY(temp_solution_tilde_algebraic, 0.5*this->timeStepWidth_, temp_increment_1);

    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    temp_solution_tilde_algebraic, temp_increment_2, timeStepNo + 1, (currentTime + 0.5*this->timeStepWidth_));

    VecAXPY(temp_increment_2, 1.0, temp_increment_1);

    VecAXPY(temp_solution_tilde, 0.25*this->timeStepWidth_, temp_increment_2);

    // second step
    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    temp_solution_tilde, temp_increment_1, timeStepNo, (currentTime + 0.5*this->timeStepWidth_));

    VecCopy(temp_solution_tilde, temp_solution_tilde_algebraic);

    VecAXPY(temp_solution_tilde_algebraic, 0.5*this->timeStepWidth_, temp_increment_1);

    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    temp_solution_tilde_algebraic, temp_increment_2, timeStepNo + 1, (currentTime + this->timeStepWidth_));

    VecAXPY(temp_increment_2, 1.0, temp_increment_1);

    VecAXPY(temp_solution_tilde, 0.25*this->timeStepWidth_, temp_increment_2);

    // check, if solutions are equal
    PetscBool flag = PETSC_FALSE;
    //VecEqual(temp_solution_tilde,temp_solution_normal, &flag);

    // calculate estimator based on current values
    if (!flag)
    {
      VecAXPY(temp_solution_tilde, -1.0, temp_solution_normal);
      VecNorm(temp_solution_tilde, NORM_2, &vecnorm_);

      PetscReal vecnorm_Combined = 0;
      MPI_Allreduce(&vecnorm_, &vecnorm_Combined, 1, MPIU_REAL, MPIU_MAX, this->data_->solution()->functionSpace()->meshPartition()->rankSubset()->mpiCommunicator());

      LOG(DEBUG) << "own vecnorm: " << vecnorm_;
      LOG(DEBUG) << ", sum of vecnorms: " << vecnorm_Combined;

      // prevent that the estimator gets nan, if the divisor is 0
      if (fabs(this->timeStepWidth_) < 1e-12)
      {
        estimator_ = 1e-10;
        LOG(DEBUG) << "estimator would be nan, set to 1e-10.";
      }
      else if (fabs(vecnorm_Combined) > 1e+75 || !std::isfinite(vecnorm_Combined))
      {
        LOG(DEBUG) << "vecnorm_Combined (" << vecnorm_Combined << ") is inf, set estimator to 1.2*tolerance_=" << tolerance_*1.2;
        estimator_ = tolerance_ * 1.2;
      }
      else
      {
        estimator_ = (double)vecnorm_ / ((1.0 - pow(0.5, 2))*this->timeStepWidth_);
        LOG(DEBUG) << "compute estimator = " << vecnorm_ << "/" << ((1.0 - pow(0.5, 2))*this->timeStepWidth_) << " = " << estimator_;
      }
    }
    else
    {
      // solutions are equal. Set low estimator to make progress and raise timeStepWidth
      estimator_ = 1e-10;
      LOG(DEBUG) << "Vecs are equal!";
    }

    // if timeStepWidth to small
    if (estimator_ <= tolerance_)
    {
      // lower multiplicator, so timeSpan is divided in lesser parts. No adjustment when multiplicator is already 1
      if (multiplicator_ > 1)
      {
        multiplicator_--;
      }

      // calculate new timeStepWidth
      this->timeStepWidth_ = timeSpan/multiplicator_;
    }
    else
    {
      // raise multiplicator, so timeSpan is divided in more parts
      if (multiplicator_ < lowestMultiplier_)
      {
        multiplicator_++;
      }

      // calculate new timeStepWidth
      this->timeStepWidth_ = timeSpan/multiplicator_;
    }

    // loop over timeSpan
    for(double time = 0.0; time < timeSpan;)
    {
      if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
      {
        // log info of current progress in timeSpan
        LOG(INFO) << "HeunAdaptive (modified), time " << time << "/" << timeSpan << ", current timeStepWidth: " << this->timeStepWidth_ << ", time remaining: " << timeSpan - time;
      }

      // calculate next solution
      this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
      solution, increment, timeStepNo, currentTime);

      VecAXPY(solution, this->timeStepWidth_, increment);

      this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
      solution, algebraicIncrement, timeStepNo + 1, currentTime + this->timeStepWidth_);

      VecAXPY(algebraicIncrement, -1.0, increment);

      VecAXPY(solution, 0.5*this->timeStepWidth_, algebraicIncrement);

      // apply the prescribed boundary condition values
      this->applyBoundaryConditions();

      // log info
      VLOG(1) << *this->data_->solution();

      // advance simulation time and timeStepNo for the Logging
      timeStepNo++;
      time = time + this->timeStepWidth_;

      // marker for the output reader
      //LOG(DEBUG) << "XXX" << currentTime << ";" << this->timeStepWidth_;

      // update loop variable and log info
      currentTime = this->startTime_ + time;

      // stop duration measurement
      if (this->durationLogKey_ != "")
        Control::PerformanceMeasurement::stop(this->durationLogKey_);

      if (!timeStepWidthsLogFilename_.empty())
      {
        // log timestep width
        timeStepWidths_.push_back(this->timeStepWidth_);
        times_.push_back(currentTime);
      }

      // write current output values
      if (withOutputWritersEnabled)
        this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);

      // start duration measurement
      if (this->durationLogKey_ != "")
        Control::PerformanceMeasurement::start(this->durationLogKey_);
    }
  }

  // log end of loop
  LOG(DEBUG) << "Exiting Loop";

  // log step info
  if (timeStepAdaptOption_ == "modified")
  {
    LOG(DEBUG) << "Step info. Successful steps: " << stepCounterSuccess_ << ", Failed steps: " << stepCounterFail_;
  }

  this->data_->solution()->restoreValuesContiguous();
  this->data_->increment()->restoreValuesContiguous();
  dataHeun->algebraicIncrement()->restoreValuesContiguous();

  // clean up to prevent memory leaks
  VecDestroy(&temp_solution_normal);
  VecDestroy(&temp_solution_tilde);
  VecDestroy(&temp_solution_tilde_algebraic);
  VecDestroy(&temp_increment_1);
  VecDestroy(&temp_increment_2);

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename DiscretizableInTime>
void HeunAdaptive<DiscretizableInTime>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTime>::run();
}

template<typename DiscretizableInTime>
double HeunAdaptive<DiscretizableInTime>::currentHeunTime()
{
  return currentTimeHeun_;
}

} // namespace TimeSteppingScheme

#include "time_stepping_scheme/heun_adaptiv.h"

#include <Python.h>
#include <memory>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTime>
HeunAdaptiv<DiscretizableInTime>::HeunAdaptiv(DihuContext context) :
  TimeSteppingExplicit<DiscretizableInTime>(context, "HeunAdaptiv")
{
  // Note, the data_ object has to be created in the initialize function

  //read allowed tolerance from settings
  tolerance_ = this->specificSettings_.getOptionDouble("tolerance", 0.1, PythonUtility::Positive);

  //read delta from settings
  delta_ = this->specificSettings_.getOptionDouble("delta", 0.1, PythonUtility::Positive);

  // backup time step width, initialized as -1 and less then timespan/2
  savedTimeStepWidth_ = -1;
}

template<typename DiscretizableInTime>
void HeunAdaptiv<DiscretizableInTime>::initialize()
{
  LOG(TRACE) << "HeunAdaptiv::initialize";

  // create data object for heun
  this->data_ = std::make_shared<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>>(this->context_);

  // call the parent initialize
  TimeSteppingSchemeOde<DiscretizableInTime>::initialize();
}

template<typename DiscretizableInTime>
void HeunAdaptiv<DiscretizableInTime>::advanceTimeSpan()
{

  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timeSpan of current step
  double timeSpan = this->endTime_ - this->startTime_;

  // Create local variables
  double alpha = 0.0;
  double estimator = 0.0;
  PetscReal vecnorm;

  // Log info
  LOG(DEBUG) << "HeunAdaptiv::advanceTimeSpan, timeSpan=" << timeSpan<< ", initial timeStepWidth=" << this->timeStepWidth_;

  // we need to cast the pointer type to the derived class. Otherwise the additional intermediateIncrement()-method of the class TimeSteppingHeun won't be there:
  std::shared_ptr<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>> dataHeun
    = std::static_pointer_cast<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>>(this->data_);

  // get vectors of all components in struct-of-array order, as needed by CellML (i.e. one long vector with [state0 state0 state0 ... state1 state1...]
  Vec &solution = this->data_->solution()->getValuesContiguous();
  Vec &increment = this->data_->increment()->getValuesContiguous();
  Vec &intermediateIncrement = dataHeun->intermediateIncrement()->getValuesContiguous();

  //Create temporary vectors for calculation
  Vec temp_solution_normal;
  Vec temp_solution_tilde;
  Vec temp_solution_tilde_intermediate;
  Vec temp_increment_1;
  Vec temp_increment_1_2;
  Vec temp_increment_2;

  // Duplicate current solutions to create same-structured vectors
  VecDuplicate(solution, &temp_solution_normal);
  VecDuplicate(solution, &temp_solution_tilde);
  VecDuplicate(solution, &temp_solution_tilde_intermediate);
  VecDuplicate(increment, &temp_increment_1);
  VecDuplicate(increment, &temp_increment_1_2);
  VecDuplicate(intermediateIncrement, &temp_increment_2);

  // Fill temporal vectors with current solution values
  VecCopy(solution, temp_solution_normal);
  VecCopy(solution, temp_solution_tilde);
  VecCopy(solution, temp_solution_tilde_intermediate);
  VecCopy(increment, temp_increment_1);
  VecCopy(increment, temp_increment_1_2);
  VecCopy(intermediateIncrement, temp_increment_2);

  // Check for savedTimeStepWidth
  if (savedTimeStepWidth_ > 0)
  {
    this->timeStepWidth_ = std::min(savedTimeStepWidth_, timeSpan*0.5);
    if (this->timeStepWidth_ == savedTimeStepWidth_){
        LOG(DEBUG) << "SavedTimeStepWidth loaded- Set timeStepWidth to:  " << this->timeStepWidth_;
    } else {
        LOG(DEBUG) << "SavedTimeStepWidth to big. Set timeStepWidth to:  " << this->timeStepWidth_;
    }

  }

  // counting time steps in one loop
  int timeStepNo = 0;

  // calculate current time
  double currentTime = this->startTime_;
  LOG(DEBUG) << "New timeSpan started. Current time: " << currentTime;

  // loop over time steps
  for(double time = 0.0; time < timeSpan;)
  {

    LOG(DEBUG) << "timeSpan progress: " << time << "/"<<timeSpan<<", timeStepWidth: " << this->timeStepWidth_ << ", time remaining: " << timeSpan - time;

    //Loop till next solution is found
    while(true){

        // Copy current solution in temporal vectors
        VecCopy(solution, temp_solution_normal);
        VecCopy(solution, temp_solution_tilde);
        VecCopy(solution, temp_solution_tilde_intermediate);

        VLOG(1) << "starting from solution: " << this->data_->solution();

        //---- Calculate next solution like always and store in temporary vector
        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        temp_solution_normal, temp_increment_1, timeStepNo, currentTime);

        VecAXPY(temp_solution_normal, this->timeStepWidth_, temp_increment_1);

        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        temp_solution_normal, temp_increment_2, timeStepNo + 1, (currentTime + this->timeStepWidth_));

        VecAXPY(temp_increment_2, -1.0, temp_increment_1);

        VecAXPY(temp_solution_normal, 0.5*this->timeStepWidth_, temp_increment_2);

        //---- Now calculate x_tilde. Use previous values as good as possible
        VecAXPY(temp_increment_2, 1.0, temp_increment_1);

        VecAXPY(temp_solution_tilde_intermediate, 0.5*this->timeStepWidth_, temp_increment_1);

        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        temp_solution_tilde_intermediate, temp_increment_1_2, timeStepNo+1, (currentTime + 0.5*this->timeStepWidth_));

        VecAXPY(temp_increment_1, 1.0, temp_increment_1_2);
        VecAXPY(temp_increment_2, 1.0, temp_increment_1_2);

        VecScale(temp_increment_1, 0.5);

        VecScale(temp_increment_2, 0.5);

        VecAXPY(temp_solution_tilde, 0.5*this->timeStepWidth_, temp_increment_1);

        VecAXPY(temp_solution_tilde, 0.5*this->timeStepWidth_, temp_increment_2);

        //Calculate estimator based on current values
        PetscBool flag;
        VecEqual(temp_solution_tilde,temp_solution_normal, &flag);
        if (!flag){
            VecAXPY(temp_solution_tilde, -1.0, temp_solution_normal);
            VecNorm(temp_solution_tilde, NORM_2, &vecnorm);
            estimator = (double)vecnorm / ((1 - pow(0.5, 2))*this->timeStepWidth_);
        }else{
            estimator = this->timeStepWidth_;
            LOG(DEBUG) << "Vecs are equal!";
            LOG(ERROR) << "Vecs are equal!";
        }

        //Calculate alpha
        alpha = 0.9 * pow((tolerance_/estimator), (1.0/3.0));
        LOG(DEBUG) << "With timeStepWidth: " << this->timeStepWidth_ << " calculated alpha: " << alpha << " and estimator: " << estimator;

        //Fallunterscheidung
        if (estimator <= tolerance_){

            // Take already calculated value with normal timestep
            VecCopy(temp_solution_normal, solution);

            // apply the prescribed boundary condition values
            this->applyBoundaryConditions();

            VLOG(1) << *this->data_->solution();

            // advance simulation time
            timeStepNo++;
            time = time + this->timeStepWidth_;

            // Adjust timeStepWidth based on alpha
            this->timeStepWidth_ = alpha*this->timeStepWidth_;
            //this->timeStepWidth_ = 1.5*this->timeStepWidth_;

            currentTime = this->startTime_ + time;
            LOG(DEBUG) << "Solution accepted. Enlarge new timeStepWidth to: " << this->timeStepWidth_;

            if (((timeSpan - time) + delta_) < this->timeStepWidth_)
            {
                if ((timeSpan - time) > 0)
                {
                    this->savedTimeStepWidth_ = this->timeStepWidth_;
                    this->timeStepWidth_ = timeSpan - time + delta_;
                    LOG(DEBUG) << "Remaining time to small for new timeStepWidth. Setting timestepwidth to " << this->timeStepWidth_;
                }
            }
            break;

        } else {
            // Reject solution and repeat loop with adjusted timeStepWidth
            this->timeStepWidth_ = alpha*this->timeStepWidth_;
            //this->timeStepWidth_ = 0.5*this->timeStepWidth_;
            LOG(DEBUG) << "Solution rejected. Retry with smaller timeStepWidth " << this->timeStepWidth_;
        }
    }

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values
    this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
    //this->data_->print();
  }
  LOG(DEBUG) << "End of timeSpan reached";

  this->data_->solution()->restoreValuesContiguous();
  this->data_->increment()->restoreValuesContiguous();
  dataHeun->intermediateIncrement()->restoreValuesContiguous();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename DiscretizableInTime>
void HeunAdaptiv<DiscretizableInTime>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTime>::run();
}
} // namespace TimeSteppingScheme

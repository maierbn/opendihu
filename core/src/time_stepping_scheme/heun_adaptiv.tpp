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
  this->data_ = std::make_shared<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>>(context);  // create data object for heun
}

template<typename DiscretizableInTime>
void HeunAdaptiv<DiscretizableInTime>::advanceTimeSpan()
{

  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  //--------------------------------------------------------------------------------------------------------------------
  // Create local variables
  double alpha = 0.0;
  double toleranz = 10;
  double estimator = 0.0;
  PetscReal vecnorm;
  //--------------------------------------------------------------------------------------------------------------------

  LOG(DEBUG) << "Heun::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // we need to cast the pointer type to the derived class. Otherwise the additional intermediateIncrement()-method of the class TimeSteppingHeun won't be there:
  std::shared_ptr<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>> dataHeun
    = std::static_pointer_cast<Data::TimeSteppingHeun<typename DiscretizableInTime::FunctionSpace, DiscretizableInTime::nComponents()>>(this->data_);

  // get vectors of all components in struct-of-array order, as needed by CellML (i.e. one long vector with [state0 state0 state0 ... state1 state1...]
  Vec &solution = this->data_->solution()->getValuesContiguous();
  Vec &increment = this->data_->increment()->getValuesContiguous();
  Vec &intermediateIncrement = dataHeun->intermediateIncrement()->getValuesContiguous();

  //--------------------------------------------------------------------------------------------------------------------
    //Create temporary vectors
    Vec temp_solution_normal = solution;
    Vec temp_solution_tilde = solution;
    Vec temp_increment_1 = increment;
    Vec temp_increment_2 = intermediateIncrement;
  //--------------------------------------------------------------------------------------------------------------------
  // loop over time steps
  double currentTime = this->startTime_;
  //for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  for(double time = 0.0; time < timeSpan;)
  {
    int timeStepNo = 0;
    //if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    //{
    //  LOG(INFO) << "Heun adaptiv, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    //}
    VLOG(1) << "starting from solution: " << this->data_->solution();

    //------------------------------------------------------------------------------------------------------------------
    //----CALCULATION OF THE TEMPORARY SOLUTION

    //---- CALCULATION FOR X_NORMAL
    //---- Create temporary vectors
    temp_solution_normal = solution;
    temp_solution_tilde = solution;

    //---- Calculate x_{i+1} like always
    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    temp_solution_normal, temp_increment_1, timeStepNo, currentTime);

    VecAXPY(temp_solution_normal, this->timeStepWidth_, temp_increment_1);

    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    temp_solution_normal, temp_increment_2, timeStepNo + 1, currentTime + this->timeStepWidth_);

    VecAXPY(temp_increment_2, -1.0, temp_increment_1);

    VecAXPY(temp_solution_normal, 0.5*this->timeStepWidth_, temp_increment_2);

    //---- NOW CALCULATE X_TILDE, atm easy way with euler explicit
    VecAXPY(temp_solution_tilde, 0.5*this->timeStepWidth_, temp_increment_1);

    this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    temp_solution_tilde, temp_increment_2, timeStepNo + 10, currentTime + this->timeStepWidth_);

    VecAXPY(temp_solution_tilde, 0.5*this->timeStepWidth_, temp_increment_2);

    //Calculate estimator
    PetscBool flag;
    VecEqual(temp_solution_tilde,temp_solution_normal, &flag);
    if (!flag){
        VecAXPY(temp_solution_tilde, -1.0, temp_solution_normal);
        VecNorm(temp_solution_tilde, NORM_2, &vecnorm);
        estimator = (double)vecnorm / ((1 - pow(0.5, 2))*this->timeStepWidth_);
    } else{
        estimator = this->timeStepWidth_;
    }
    //Calculate alpha
    //std::cout<<"Estimator="<<estimator<<"\n";
    alpha = pow((toleranz/estimator), (1/3));

    //Fallunterscheidung
    if (estimator <= toleranz){
        // advance solution value to compute u* first
        // compute  delta_u = f(u_{t})
        // we call f(u_{t}) the "increment"
        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        solution, increment, timeStepNo, currentTime);

        // integrate u* += dt * delta_u : values = solution.values + timeStepWidth * increment.values
        VecAXPY(solution, this->timeStepWidth_, increment);

        VLOG(1) << "increment: " << this->data_->increment() << ", dt: " << this->timeStepWidth_;

        // now, advance solution value to compute u_{t+1}
        // compute  delta_u* = f(u*)
        // we call f(u*) the "intermediateIncrement"
        this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
        solution, intermediateIncrement, timeStepNo + 1, currentTime + this->timeStepWidth_);

        // integrate u_{t+1} = u_{t} + dt*0.5(delta_u + delta_u_star)
        // however, use: u_{t+1} = u* + 0.5*dt*(f(u*)-f(u_{t}))     (#)
        //
        // first calculate (f(u*)-f(u_{t})). to save storage we store into f(u*):
        VecAXPY(intermediateIncrement, -1.0, increment);

        // now compute overall step as described above (#)
        VecAXPY(solution, 0.5*this->timeStepWidth_, intermediateIncrement);

        // apply the prescribed boundary condition values
        this->applyBoundaryConditions();

        VLOG(1) << *this->data_->solution();

        // advance simulation time
        //timeStepNo++;
        time = time + this->timeStepWidth_;
        //std::cout<<"Time = " << time;
        //time = time + timeSpan;
        this->timeStepWidth_ = this->timeStepWidth_*alpha;
        //std::cout<<"Timestep accepted, alpha = " << alpha << " , next timestep: "<<this->timeStepWidth_<<"\n";
        currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    } else {
        this->timeStepWidth_ = this->timeStepWidth_*alpha;
        //std::cout<<"Timestep rejected, next timestep: "<<this->timeStepWidth_<<"\n";
    }

    //------------------------------------------------------------------------------------------------------------------

    // advance solution value to compute u* first
    // compute  delta_u = f(u_{t})
    // we call f(u_{t}) the "increment"
    //this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    //  solution, increment, timeStepNo, currentTime);

    // integrate u* += dt * delta_u : values = solution.values + timeStepWidth * increment.values
    //VecAXPY(solution, this->timeStepWidth_, increment);

    //VLOG(1) << "increment: " << this->data_->increment() << ", dt: " << this->timeStepWidth_;

    // now, advance solution value to compute u_{t+1}
    // compute  delta_u* = f(u*)
    // we call f(u*) the "intermediateIncrement"
    //this->discretizableInTime_.evaluateTimesteppingRightHandSideExplicit(
    //  solution, intermediateIncrement, timeStepNo + 1, currentTime + this->timeStepWidth_);

    // integrate u_{t+1} = u_{t} + dt*0.5(delta_u + delta_u_star)
    // however, use: u_{t+1} = u* + 0.5*dt*(f(u*)-f(u_{t}))     (#)
    //
    // first calculate (f(u*)-f(u_{t})). to save storage we store into f(u*):
    //VecAXPY(intermediateIncrement, -1.0, increment);

    // now compute overall step as described above (#)
    //VecAXPY(solution, 0.5*this->timeStepWidth_, intermediateIncrement);

    // apply the prescribed boundary condition values
    //this->applyBoundaryConditions();

    //VLOG(1) << *this->data_->solution();

    // advance simulation time
    //timeStepNo++;
    //time = time + timeSpan;
    //currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

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

#include "time_stepping_scheme/crank_nicolson.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include <petscksp.h>
#include "solver/solver_manager.h"
#include "solver/linear.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
CrankNicolson<DiscretizableInTimeType>::CrankNicolson(DihuContext context) :
TimeSteppingImplicit<DiscretizableInTimeType>(context, "CrankNicolson")
{
}

template<typename DiscretizableInTimeType>
void CrankNicolson<DiscretizableInTimeType>::advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;
  
  LOG(DEBUG) << "CrankNicolson::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  Vec solution = this->data_->solution()->valuesGlobal();
  Vec systemRightHandSide = this->dataImplicit_->systemRightHandSide()->valuesGlobal();

  // loop over time steps
  double currentTime = this->startTime_;
  
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "Crank Nicolson, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }
    
    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    
    // compute systemRightHandSide = solution * integrationMatrix
    this->evaluateTimesteppingRightHandSideImplicit(solution, systemRightHandSide, timeStepNo, currentTime);

    // adjust rhs vector such that boundary conditions are satisfied
    this->dirichletBoundaryConditions_->applyInRightHandSide(this->dataImplicit_->systemRightHandSide(), this->dataImplicit_->boundaryConditionsRightHandSideSummand());

    // advance computed value
    // solve A*u^{t+1} = u^{t} for u^{t+1} where A is the system matrix, solveLinearSystem(b,x)
    this->solveLinearSystem(systemRightHandSide, solution);
    
    // check if the solution contains Nans or Inf values
    this->checkForNanInf(timeStepNo, currentTime);

    VLOG(1) << *this->data_->solution();

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

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename DiscretizableInTimeType>
void CrankNicolson<DiscretizableInTimeType>::
initialize()
{
  LOG(TRACE) << "CrankNicolson::initialize()";
  
  if (this->initialized_)
    return;
  
  TimeSteppingImplicit<DiscretizableInTimeType>::initialize();
  LOG(TRACE) << "setIntegrationMatrixRightHandSide()";
  this->setIntegrationMatrixRightHandSide();
  
  this->initialized_ = true;
}

template<typename DiscretizableInTimeType>
void CrankNicolson<DiscretizableInTimeType>::
setSystemMatrix(double timeStepWidth)
{
  LOG(TRACE) << "setSystemMatrix(timeStepWidth=" << timeStepWidth << ")";

  // compute the system matrix (I - dt*M^{-1}K) where M^{-1} is the lumped mass matrix
  
  Mat &inverseLumpedMassMatrix = this->discretizableInTime_.data().inverseLumpedMassMatrix()->valuesGlobal();
  Mat &stiffnessMatrix = this->discretizableInTime_.data().stiffnessMatrix()->valuesGlobal();
  Mat systemMatrix;
  
  PetscErrorCode ierr;
  
  // compute systemMatrix = M^{-1}K
  // the result matrix is created by MatMatMult
  ierr = MatMatMult(inverseLumpedMassMatrix, stiffnessMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &systemMatrix);
  this->dataImplicit_->initializeSystemMatrix(systemMatrix);
  
  // scale systemMatrix by -dt, systemMatrix = -dt/2 *M^{-1}K
  ierr = MatScale(this->dataImplicit_->systemMatrix()->valuesGlobal(), -0.5*timeStepWidth); CHKERRV(ierr);
  
  // add 1 on the diagonal: systemMatrix = I - dt/2 *M^{-1}K
  ierr = MatShift(this->dataImplicit_->systemMatrix()->valuesGlobal(), 1.0); CHKERRV(ierr);
  
  this->dataImplicit_->systemMatrix()->assembly(MAT_FINAL_ASSEMBLY);
  
  VLOG(1) << *this->dataImplicit_->systemMatrix();
}

template<typename DiscretizableInTimeType>
void CrankNicolson<DiscretizableInTimeType>::
setIntegrationMatrixRightHandSide()
{
  LOG(TRACE) << "setIntegrationMatrixRightHandSide()";

  Mat &systemMatrix = this->dataImplicit_->systemMatrix()->valuesGlobal();
  Mat integrationMatrix; 
  
  PetscErrorCode ierr;
  
  // copy integration matrix from the system matrix
  ierr = MatConvert(systemMatrix, MATSAME, MAT_INITIAL_MATRIX, &integrationMatrix); CHKERRV(ierr); //it creates the new matrix
  
  // scale systemMatrix by -dt, systemMatrix = dt/2*M^{-1}K
  ierr = MatScale(integrationMatrix, -1.0); CHKERRV(ierr);
  
  // add 1 on the diagonal: systemMatrix = I + dt/2*M^{-1}K
  ierr = MatShift(integrationMatrix, 2.0); CHKERRV(ierr);
  
  this->dataImplicit_->initializeIntegrationMatrixRightHandSide(integrationMatrix);
  
  //this->dataImplicit_->initializeMatrix(integrationMatrix, this->dataImplicit_->integrationMatrixRightHandSide(), "integrationMatrixRightHandSide");
  
  //qierr=MatView(this->dataImplicit_->integrationMatrixRightHandSide()->valuesGlobal(), PETSC_VIEWER_STDOUT_WORLD);
  
  this->dataImplicit_->integrationMatrixRightHandSide()->assembly(MAT_FINAL_ASSEMBLY);
  
  VLOG(1) << *this->dataImplicit_->integrationMatrixRightHandSide();

}

template<typename DiscretizableInTimeType>
void CrankNicolson<DiscretizableInTimeType>::
evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  // LOG(TRACE) << "evaluateTimesteppingRightHandSideImplicit";
  
  // this method computes output = input * (I+dt/2 M^(-1) K)= input *(-A+2I), where A=(I-dt/2 M^(-1) K) is the system matrix
  Mat &integrationMatrix = this->dataImplicit_->integrationMatrixRightHandSide()->valuesGlobal();
  
  PetscErrorCode ierr;
  ierr = MatMult(integrationMatrix, input, output); CHKERRV(ierr);    // MatMult(mat,x,y) computes y = Ax
  
  VLOG(1) << PetscUtility::getStringVector(output);
}

} // namespace TimeSteppingScheme

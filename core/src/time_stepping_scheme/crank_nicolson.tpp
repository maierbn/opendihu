#include "time_stepping_scheme/implicit_euler.h"

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
TimeSteppingImplicit<DiscretizableInTimeType>(context)
{
  //this->data_ = std::make_shared <Data::TimeStepping<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>>(context); // create data object for implicit scheme
  PyObject *topLevelSettings = this->context_.getPythonConfig();
  this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "CrankNicolson");
  this->outputWriterManager_.initialize(this->specificSettings_);
}

template<typename DiscretizableInTimeType>
void CrankNicolson<DiscretizableInTimeType>::advanceTimeSpan()
{
  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;
  
  LOG(DEBUG) << "CrankNicolson::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
  << " n steps: " << this->numberTimeSteps_;
  
  int nEntries;
  VecGetSize(this->data_->solution().valuesGlobal(), &nEntries);
  
  // loop over time steps
  double currentTime = this->startTime_;
  
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0)
    {
      std::stringstream threadNumberMessage;
      threadNumberMessage << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
      LOG(INFO) << threadNumberMessage.str() << ": Timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }
    
    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
    
    Vec rhs;
    this->evaluateTimesteppingRightHandSideImplicit(this->data_->solution().valuesGlobal(), &rhs, timeStepNo, currentTime)
    
    // advance computed value
    // solve A*u^{t+1} = u^{t} for u^{t+1} where A is the system matrix
    
    this->solveLinearSystem(&rhs, this->data_->solution().valuesGlobal());
    
    // write current output values
    this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime);
    
    //this->data_->print();
  } 
}

template<typename DiscretizableInTimeType>
void CrankNicolson<DiscretizableInTimeType>::
initialize()
{
  TimeSteppingImplicit<DiscretizableInTimeType>::initialize();
  
  this->setIntegrationMatrixRightHandSide();
  
  this->initialized_ = true;
}

template<typename DiscretizableInTimeType>
void CrankNicolson<DiscretizableInTimeType>::
setSystemMatrix(double timeStepWidth)
{
  LOG(TRACE) << "setSystemMatrix(timeStepWidth=" << timeStepWidth << ")";
  
  //if(!this->discretizableInTime_.invLumMassMatrixSet())
    //this->discretizableInTime_.setInverseLumpedMassMatrix();
  
  // compute the system matrix (I - dt*M^{-1}K) where M^{-1} is the lumped mass matrix
  
  Mat &inverseLumpedMassMatrix = this->discretizableInTime_->data_.inverseLumpedMassMatrix()->valuesGlobal();
  Mat &stiffnessMatrix = this->discretizableInTime_->data_.stiffnessMatrix()->valuesGlobal();
  Mat systemMatrix;
  
  PetscErrorCode ierr;
  
  // compute systemMatrix = M^{-1}K
  // the result matrix is created by MatMatMult
  ierr = MatMatMult(inverseLumpedMassMatrix, stiffnessMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &systemMatrix);
  this->data_.initializeSystemMatrix(systemMatrix);
  
  // scale systemMatrix by -dt, systemMatrix = -dt/2 *M^{-1}K
  ierr = MatScale(this->data_.systemMatrix()->valuesGlobal(), -0.5*timeStepWidth); CHKERRV(ierr);
  
  // add 1 on the diagonal: systemMatrix = I - dt/2 *M^{-1}K
  ierr = MatShift(this->data_.systemMatrix()->valuesGlobal(), 1.0); CHKERRV(ierr);
  
  this->data_.systemMatrix()->assembly(MAT_FINAL_ASSEMBLY);
  
  VLOG(1) << *this->data_.systemMatrix();
}

template<typename DiscretizableInTimeType>
void CrankNicolson<DiscretizableInTimeType>::
setIntegrationMatrixRightHandSide()
{
  LOG(TRACE) << "setIntegrationMatrixRightHandSideImplicit()";
  
  Mat &systemMatrix = this->data_.systemMatrix()->valuesGlobal();
  Mat integrationMatrix= this->data_.integrationMatrixRightHandSide()->valuesGlobal();
  
  PetscErrorCode ierr;
  
  // copy integration matrix from the system matrix
  ierr=MatDuplicate(systemMatrix, MAT_COPY_VALUES, integrationMatrix); CHKERRV(ierr); //does it initializes the object? 
  
  // scale systemMatrix by -dt, systemMatrix = dt/2*M^{-1}K
  ierr = MatScale(this->data_.integrationMatrixRightHandSide()->valuesGlobal(), -1.0); CHKERRV(ierr);
  
  // add 1 on the diagonal: systemMatrix = I + dt/2*M^{-1}K
  ierr = MatShift(this->data_.integrationMatrixRightHandSide()->valuesGlobal(), 2.0); CHKERRV(ierr);
  
  this->data_.integrationMatrixRightHandSide()->assembly(MAT_FINAL_ASSEMBLY);
  
  VLOG(1) << *this->data_.integrationMatrixRightHandSide();
}

template<DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  LOG(TRACE) << "evaluateTimesteppingRightHandSideImplicit";
  
  // this method computes output = input * (I+dt/2 M^(-1) K)= input *(-A+2I), where A=(I-dt/2 M^(-1) K) is the system matrix
  
  Mat &integrationMatrix= this->data_.integrationMatrixRightHandSide()->valuesGlobal();
  
  PetscErrorCode ierr;
  
  ierr = MatMult(integrationMatrix, input, output); CHKERRV(ierr);
  
  VLOG(1) << PetscUtility::getStringVector(output);
}

} // namespace TimeSteppingScheme

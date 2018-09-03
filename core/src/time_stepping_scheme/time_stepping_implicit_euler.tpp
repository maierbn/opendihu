#include "time_stepping_scheme/time_stepping_implicit_euler.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include <petscksp.h>
#include "solver/solver_manager.h"
#include "solver/linear.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
ImplicitEuler<DiscretizableInTimeType>::ImplicitEuler(DihuContext context) :
TimeSteppingImplicit<DiscretizableInTimeType>(context)
{
  //this->data_ = std::make_shared <Data::TimeStepping<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>>(context); // create data object for implicit euler
  //PyObject *topLevelSettings = this->context_.getPythonConfig();
  //this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "ImplicitEuler");
  //this->outputWriterManager_.initialize(this->specificSettings_);
}

/*
template<typename DiscretizableInTimeType>
void ImplicitEuler<DiscretizableInTimeType>::advanceTimeSpan()
{
  TimeSteppingImplicit<DiscretizableInTimeType>::advanceTimeSpan();  
}

template<typename DiscretizableInTimeType>
void ImplicitEuler<DiscretizableInTimeType>::run()
{
  TimeSteppingSchemeOde<DiscretizableInTimeType>::run();
}
*/

template<typename DiscretizableInTimeType>
void ImplicitEuler<DiscretizableInTimeType>::
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
  
  // scale systemMatrix by -dt, systemMatrix = -dt*M^{-1}K
  ierr = MatScale(this->data_.systemMatrix()->valuesGlobal(), -timeStepWidth); CHKERRV(ierr);
  
  // add 1 on the diagonal: systemMatrix = I - dt*M^{-1}K
  ierr = MatShift(this->data_.systemMatrix()->valuesGlobal(), 1.0); CHKERRV(ierr);
  
  this->data_.systemMatrix()->assembly(MAT_FINAL_ASSEMBLY);
  
  VLOG(1) << *this->data_.systemMatrix();
}

} // namespace TimeSteppingScheme

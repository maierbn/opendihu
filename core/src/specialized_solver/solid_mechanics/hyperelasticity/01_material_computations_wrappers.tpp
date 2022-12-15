#include "specialized_solver/solid_mechanics/hyperelasticity/01_material_computations.h"

#include <Python.h>  // has to be the first included header

namespace SpatialDiscretization
{

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
bool HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
materialComputeInternalVirtualWork(
  std::shared_ptr<VecHyperelasticity> displacements,
  std::shared_ptr<VecHyperelasticity> internalVirtualWork
)
{
  // this is only for the static problem
  assert(nDisplacementComponents == 3);
  assert(displacements);
  assert(internalVirtualWork);

  // copy the values from the combined vector to displacements and pressure, this is needed for materialComputeResidual
  setDisplacementsAndPressureFromCombinedVec(displacements->valuesGlobal());

  // place the result vector to be used by materialComputeInternalVirtualWork
  Vec result = internalVirtualWork->valuesGlobal();
  PetscErrorCode ierr;

  ierr = VecSwap(result, solverVariableResidual_);
  if(ierr)
    return false;

  // compute the actual output of the nonlinear function
  bool successful = materialComputeInternalVirtualWork();

  // reset the result vector
  ierr = VecSwap(result, solverVariableResidual_);
  if(ierr)
    return false;
  return successful;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
solveForDisplacements(
  std::shared_ptr<VecHyperelasticity> externalVirtualWork,
  std::shared_ptr<VecHyperelasticity> displacements
)
{
  // solve ∂W_int - ∂W_ext = 0 with J = 1

  // assign new value for ∂W_ext, save previous vector
  Vec rhs = externalVirtualWork->valuesGlobal();
  PetscErrorCode ierr;
  ierr = VecSwap(rhs, externalVirtualWorkDead_); CHKERRV(ierr);

  Vec result = displacements->valuesGlobal();
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);

  this->nonlinearSolve();

  // restore value of ∂W_ext and solution
  ierr = VecSwap(rhs, externalVirtualWorkDead_); CHKERRV(ierr);
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
solveDynamicProblem(
  std::shared_ptr<VecHyperelasticity> displacementsVelocitiesPressure, bool isFirstTimeStep,
  Vec internalVirtualWork, Vec &externalVirtualWorkDead, Vec accelerationTerm, bool withOutputWritersEnabled
)
{
  // solve δW_int - δW_ext,dead + ∫_Ω ρ0 1/dt (v^(n+1),L - v^(n),L)/dt ϕ^L ϕ^M δu dV = 0,
  // 1/dt (u^(n+1) - u^(n)) - v^(n+1) = 0
  // ∫ (J-1) δp = 0   ∀ δp,δu,δv

  assert(nDisplacementComponents == 6);

  // use given solution vector
  PetscErrorCode ierr;
  Vec result = displacementsVelocitiesPressure->valuesGlobal();
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);

  // write reference output values
  if (isFirstTimeStep)
  {
    // write initial values but don't increment counter
    this->outputWriterManager_.writeOutput(this->data_, 0, 0.0, 0);
    this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, 0.0, 0);
  }

  if (this->extrapolateInitialGuess_)
  {
    // copy the solution values back to this->data_.displacements(), and this->data_.velocities() and this->data.pressure()
    setUVP(combinedVecSolution_->valuesGlobal());

    // determine new initial values by extrapolating u,v over timestep width
    Vec u = this->data_.displacements()->valuesGlobal();
    Vec v = this->data_.velocities()->valuesGlobal();
    Vec vprev = this->data_.velocitiesPreviousTimestep()->valuesGlobal();

    LOG(DEBUG) << "extrapolateInitialGuess: initial u: " << *this->data_.displacements();
    LOG(DEBUG) << "extrapolateInitialGuess: initial v: " << *this->data_.velocities();

    // displacements = u + dt*v
    ierr = VecAXPY(u, this->timeStepWidth_, v); CHKERRV(ierr);

    // compute a = (v^(i) - v^(i-1)) / dt in the buffer velocitiesPreviousTimestep()
    // velocitiesPreviousTimestep = v^(i-1) - v^(i)
    ierr = VecAXPY(vprev, -1, v); CHKERRV(ierr);

    // velocitiesPreviousTimestep = (v^(i) - v^(i-1)) / dt = acceleration
    ierr = VecScale(vprev, -1.0/this->timeStepWidth_); CHKERRV(ierr);
    Vec a = vprev;

    // velocities = v + dt*a
    ierr = VecAXPY(v, this->timeStepWidth_, a); CHKERRV(ierr);

    LOG(DEBUG) << "extrapolateInitialGuess: new u: " << *this->data_.displacements();
    LOG(DEBUG) << "extrapolateInitialGuess: new v: " << *this->data_.velocities();
    LOG(DEBUG) << "extrapolateInitialGuess: a: " << *this->data_.velocitiesPreviousTimestep();
    LOG(DEBUG) << "dt: " << this->timeStepWidth_;
  }

  // save u,v,p of previous timestep from solverVariableSolution_
  this->setDisplacementsVelocitiesAndPressureFromCombinedVec(solverVariableSolution_,
                                                             this->data_.displacementsPreviousTimestep(),
                                                             this->data_.velocitiesPreviousTimestep(),
                                                             this->data_.pressurePreviousTimestep());

  LOG(DEBUG) << "saved initial uvp values: ";
  LOG(DEBUG) << *this->data_.displacementsPreviousTimestep();
  VLOG(1) << *this->data_.velocitiesPreviousTimestep();
  VLOG(1) << *this->data_.pressurePreviousTimestep();

  if (this->extrapolateInitialGuess_)
  {
    combinedVecSolution_->startGhostManipulation();

    static std::vector<double> values;
    // store new displacements in combinedVecSolution_ which is solverVariableSolution_
    for (int i = 0; i < 3; i++)
    {
      values.clear();
      this->data_.displacements()->getValuesWithoutGhosts(i, values);
      combinedVecSolution_->setValues(i, displacementsFunctionSpace_->meshPartition()->nDofsLocalWithoutGhosts(),
                                      displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data(), INSERT_VALUES);
    }

    // store new velocities in combinedVecSolution_
    for (int i = 0; i < 3; i++)
    {
      values.clear();
      this->data_.velocities()->getValuesWithoutGhosts(i, values);
      combinedVecSolution_->setValues(3+i, displacementsFunctionSpace_->meshPartition()->nDofsLocalWithoutGhosts(),
                                      displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data(), INSERT_VALUES);
    }

    combinedVecSolution_->zeroGhostBuffer();
    combinedVecSolution_->finishGhostManipulation();

    this->data_.displacements()->zeroGhostBuffer();
    this->data_.displacements()->finishGhostManipulation();
    this->data_.velocities()->zeroGhostBuffer();
    this->data_.velocities()->finishGhostManipulation();

    this->data_.displacements()->startGhostManipulation();
    this->data_.velocities()->startGhostManipulation();

    LOG(DEBUG) << "extrapolateInitialGuess: initial values: " << getString(solverVariableSolution_);
  }

  // find the solution for u,v,p of the nonlinear equation, potentially multiple load steps and repeated solve calls
  this->nonlinearSolve();

  LOG(DEBUG) << "result: " << getString(solverVariableSolution_);

  // update geometry fields from displacements, compute PK2Stress, write output with output writers
  this->postprocessSolution();

  LOG(DEBUG) << "nonlinearSolve finished";

  std::vector<Vec3> displacementValues;
  this->data_.displacements()->getValuesWithoutGhosts(displacementValues);

  // output with output writers of hyperelasticity_solver
  if (withOutputWritersEnabled)
  {
    this->outputWriterManager_.writeOutput(this->data_, 0, this->endTime_);
    this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, this->endTime_);
  }

  // determine output values for δW
  externalVirtualWorkDead = externalVirtualWorkDead_;

  // compute internal virtual work
  setUVP(solverVariableSolution_);
  materialComputeInternalVirtualWork();
  ierr = VecCopy(solverVariableResidual_, internalVirtualWork); CHKERRV(ierr);

  // compute acceleration term
  ierr = VecZeroEntries(solverVariableResidual_); CHKERRV(ierr);
  materialAddAccelerationTermAndVelocityEquation();
  ierr = VecCopy(solverVariableResidual_, accelerationTerm); CHKERRV(ierr);

  // restore value of ∂W_ext and solution
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityMaterialComputations<Term,withLargeOutput,MeshType,nDisplacementComponents>::
solveQuasistaticProblem(
  std::shared_ptr<VecHyperelasticity> displacementsVelocitiesPressure, bool isFirstTimeStep, bool withOutputWritersEnabled
)
{
  // solve δW_int - δW_ext,dead + ∫_Ω ρ0 1/dt (v^(n+1),L - v^(n),L)/dt ϕ^L ϕ^M δu dV = 0,
  // 1/dt (u^(n+1) - u^(n)) - v^(n+1) = 0
  // ∫ (J-1) δp = 0   ∀ δp,δu,δv

  assert(nDisplacementComponents == 6);

  // use given solution vector
  PetscErrorCode ierr;
  Vec result = displacementsVelocitiesPressure->valuesGlobal();
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);

  // write reference output values
  if (isFirstTimeStep)
  {
    // write initial values but don't increment counter
    this->outputWriterManager_.writeOutput(this->data_, 0, 0.0, 0);
    this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, 0.0, 0);
  }

  // save u,v,p of previous timestep from solverVariableSolution_
  this->setDisplacementsVelocitiesAndPressureFromCombinedVec(solverVariableSolution_,
                                                             this->data_.displacementsPreviousTimestep(),
                                                             this->data_.velocitiesPreviousTimestep(),
                                                             this->data_.pressurePreviousTimestep());

  LOG(DEBUG) << "saved initial uvp values: ";
  LOG(DEBUG) << *this->data_.displacementsPreviousTimestep();
  VLOG(1) << *this->data_.velocitiesPreviousTimestep();
  VLOG(1) << *this->data_.pressurePreviousTimestep();

  // find the solution for u,v,p of the nonlinear equation, potentially multiple load steps and repeated solve calls
  this->nonlinearSolve();

  LOG(DEBUG) << "result: " << getString(solverVariableSolution_);

  // update geometry fields from displacements, compute PK2Stress, write output with output writers
  this->postprocessSolution();

  LOG(DEBUG) << "nonlinearSolve finished";

  // output with output writers of hyperelasticity_solver
  if (withOutputWritersEnabled)
  {
    this->outputWriterManager_.writeOutput(this->data_, 0, this->endTime_);
    this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, this->endTime_);
  }

}

}  // namespace SpatialDiscretization

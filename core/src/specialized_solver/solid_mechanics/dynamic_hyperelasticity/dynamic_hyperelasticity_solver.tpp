#include "specialized_solver/solid_mechanics/dynamic_hyperelasticity/dynamic_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.h"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"

namespace TimeSteppingScheme
{

template<typename Term>
DynamicHyperelasticitySolver<Term>::
DynamicHyperelasticitySolver(DihuContext context) :
  TimeSteppingScheme(context["DynamicHyperelasticitySolver"]), staticSolver_(this->context_), data_(context_)
{
  specificSettings_ = this->context_.getPythonConfig();
  density_ = specificSettings_.getOptionDouble("density", 1.0, PythonUtility::Positive);
  viscosity_ = specificSettings_.getOptionDouble("viscosity", 0.0, PythonUtility::NonNegative);
}

template<typename Term>
void DynamicHyperelasticitySolver<Term>::
initialize()
{
  TimeSteppingScheme::initialize();
  staticSolver_.initialize();
  data_.setFunctionSpace(staticSolver_.data().functionSpace());
  data_.initialize(staticSolver_.data().displacements());

  // initialize temporary vectors
  u_ = staticSolver_.createPartitionedPetscVec("u");
  v_ = staticSolver_.createPartitionedPetscVec("v");
  a_ = staticSolver_.createPartitionedPetscVec("a");

  for (int i = 0; i < 9; i++)
  {
    std::stringstream name;
    name << "temp" << i;
    temp_[i] = staticSolver_.createPartitionedPetscVec(name.str());
  }

  for (int i = 0; i < 4; i++)
  {
    std::stringstream name;
    name << "k" << i;
    k_[i] = staticSolver_.createPartitionedPetscVec(name.str());

    name.str();
    name << "l" << i;
    l_[i] = staticSolver_.createPartitionedPetscVec(name.str());
  }

  initializeMassMatrix();
}

template<typename Term>
void DynamicHyperelasticitySolver<Term>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "DynamicHyperelasticitySolver::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // get PETSc vectors to work on
  //Vec &displacements = this->data_.displacements()->valuesGlobal();
  //Vec &velocity = this->data_.velocity()->valuesGlobal();
  //Vec &acceleration = this->data_.acceleration()->valuesGlobal();

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "DynamicHyperelasticitySolver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    computeRungeKutta4();

    staticSolver_.run();

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values
    this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
    //this->data_.print();
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 1, endTime_);
}

template<typename Term>
void DynamicHyperelasticitySolver<Term>::
initializeMassMatrix()
{
  LOG(TRACE) << "init mass matrix";

  // create inverseLumpedMassMatrix
  massMatrix_ = staticSolver_.createPartitionedPetscMat("massMatrix");

  const int D = 3;  // dimension
  const int nComponents = 3;

  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,Quadrature::Gauss<3>> QuadratureDD;   // quadratic*quadratic = 4th order polynomial, 3 gauss points = 2*3-1 = 5th order exact

  const int nDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
  typedef MathUtility::Matrix<nDofsPerElement,nDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]

  // initialize variables
  std::shared_ptr<DisplacementsFunctionSpace> functionSpace = this->data_.functionSpace();
  functionSpace->geometryField().setRepresentationGlobal();
  functionSpace->geometryField().startGhostManipulation();   // ensure that local ghost values of geometry field are set

  // initialize values to zero
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < functionSpace->nElementsLocal(); elementNo++)
  {
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNo);

    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int j = 0; j < nDofsPerElement; j++)
      {
        for (int componentNo = 0; componentNo < nComponents; componentNo++)
        {
          massMatrix_->setValue(componentNo, dofNosLocal[i], componentNo, dofNosLocal[j], 0, INSERT_VALUES);
        }
      }
    }
  }
  massMatrix_->assembly(MAT_FLUSH_ASSEMBLY);

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  // set entries in massMatrix_
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < functionSpace->nElementsLocal(); elementNo++)
  {
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNo);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,DisplacementsFunctionSpace::nDofsPerElement()> geometry;
    functionSpace->getElementGeometry(elementNo, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = DisplacementsFunctionSpace::computeJacobian(geometry, xi);

      // get evaluations of integrand which is defined in another class
      evaluationsArray[samplingPointIndex] = SpatialDiscretization::IntegrandMassMatrix<D,EvaluationsType,DisplacementsFunctionSpace,1,Term>::evaluateIntegrand(jacobian,xi) * density_;

    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration and add to entry in rhs vector
    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int j = 0; j < nDofsPerElement; j++)
      {
        // integrate value and set same entry for all components
        double integratedValue = integratedValues(i, j);

        for (int componentNo = 0; componentNo < nComponents; componentNo++)
        {
          massMatrix_->setValue(componentNo, dofNosLocal[i], componentNo, dofNosLocal[j], integratedValue, ADD_VALUES);
        }
      }  // j
    }  // i
  }  // elementNo

  // merge local changes in parallel and assemble the matrix (MatAssemblyBegin, MatAssemblyEnd)
  massMatrix_->assembly(MAT_FINAL_ASSEMBLY);
  LOG(DEBUG) << *massMatrix_;

  inverseLumpedMassMatrix_ = staticSolver_.createPartitionedPetscMat("inverseLumpedMassMatrix");

  // row sum vector
  std::shared_ptr<VecHyperelasticity> rowSum = staticSolver_.createPartitionedPetscVec("rowSum");

  LOG(DEBUG) << rowSum;

  // store the sum of each row of the matrix in the vector rowSum
  PetscErrorCode ierr;
  ierr = MatGetRowSum(massMatrix_->valuesGlobal(), rowSum->valuesGlobal()); CHKERRV(ierr);

  // for the inverse matrix, replace each entry in rowSum by its reciprocal
  ierr = VecReciprocal(rowSum->valuesGlobal()); CHKERRV(ierr);

  // set the values on the diagonal
  ierr = MatDiagonalSet(inverseLumpedMassMatrix_->valuesGlobal(), rowSum->valuesGlobal(), INSERT_VALUES); CHKERRV(ierr);

  this->inverseLumpedMassMatrix_->assembly(MAT_FINAL_ASSEMBLY);

  LOG(DEBUG) << *inverseLumpedMassMatrix_;
}

template<typename Term>
void DynamicHyperelasticitySolver<Term>::
computeAcceleration(std::shared_ptr<VecHyperelasticity> u, std::shared_ptr<VecHyperelasticity> v, std::shared_ptr<VecHyperelasticity> a)
{
  assert(u);
  assert(v);
  assert(a);
  LOG(TRACE) << "computeAcceleration";

  // compute δW_int, store in temp_[0]
  staticSolver_.materialComputeInternalVirtualWork(u, temp_[0]);

  // subtract external virtual work f

  // δW_ext = int_∂Ω T_a phi_L dS was precomputed in initialize of staticSolver_
  PetscErrorCode ierr;
  ierr = VecAXPY(temp_[0]->valuesGlobal(), -1, staticSolver_.externalVirtualWork()); CHKERRV(ierr);

  LOG(DEBUG) << "-f + Ku: " << *temp_[0];

  // compute damping and add to temp, then temp0 = -f + Ku + Cv
  temp_[1]->zeroEntries();
  addDamping(v, temp_[1]);

  LOG(DEBUG) << "damping: " << *temp_[1];

  ierr = VecAXPY(temp_[0]->valuesGlobal(), 1, temp_[1]->valuesGlobal()); CHKERRV(ierr);

  LOG(DEBUG) << "-f + Ku + Cv: " << *temp_[0];

  // negate, then temp0 = f - Ku - Cv
  ierr = VecScale(temp_[0]->valuesGlobal(), -1); CHKERRV(ierr);

  // premultiply inverse lumped mass matrix, a = M^-1(temp0) = M^-1*(f - Ku - Cv)
  ierr = MatMult(inverseLumpedMassMatrix_->valuesGlobal(), temp_[0]->valuesGlobal(), a->valuesGlobal()); CHKERRV(ierr);
}

template<typename Term>
void DynamicHyperelasticitySolver<Term>::
addDamping(std::shared_ptr<VecHyperelasticity> v, std::shared_ptr<VecHyperelasticity> damping)
{
  LOG(TRACE) << "addDamping";
  // compute the damping
  // C*v with C_ij = ∫φ_i*μ*φ_j dV (for all components), v = velocity
  // damping_LbMb = ∫φ_Lb*μ*φ_Mb*v_Mb dV (no sum over b, but over M), damping_LaMb = 0 for a != b

  const int D = 3;  // dimension
  const int nComponents = 3;

  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,Quadrature::Gauss<3>> QuadratureDD;   // quadratic*quadratic = 4th order polynomial, 3 gauss points = 2*3-1 = 5th order exact

  const int nDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
  typedef MathUtility::Matrix<nDofsPerElement,nDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]

  // initialize variables
  std::shared_ptr<DisplacementsFunctionSpace> functionSpace = this->data_.functionSpace();
  functionSpace->geometryField().setRepresentationGlobal();
  functionSpace->geometryField().startGhostManipulation();   // ensure that local ghost values of geometry field are set

  damping->setRepresentationGlobal();
  damping->startGhostManipulation();
  damping->zeroGhostBuffer();

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  // set entries in massMatrix_
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < functionSpace->nElementsLocal(); elementNo++)
  {
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNo);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,DisplacementsFunctionSpace::nDofsPerElement()> geometry;
    functionSpace->getElementGeometry(elementNo, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = DisplacementsFunctionSpace::computeJacobian(geometry, xi);

      // get evaluations of integrand which is defined in another class
      evaluationsArray[samplingPointIndex]
        = SpatialDiscretization::IntegrandMassMatrix<D,EvaluationsType,DisplacementsFunctionSpace,1,Term>::evaluateIntegrand(jacobian,xi) * viscosity_;

    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration and add to entry in vector
    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int j = 0; j < nDofsPerElement; j++)
      {
        // integrate value: ∫φ_i*μ*φ_j dV
        double integratedValue = integratedValues(i, j);

        for (int componentNo = 0; componentNo < nComponents; componentNo++)
        {
          // get the velocity v_ja  (with a: componentNo)
          double velocity = v->getValue(componentNo, dofNosLocal[j]);

          // compute the damping value: damping_ia += ∫φ_i*μ* Σ_j φ_ja*v_ja dV (with i,j: dofNos and a: componentNo)
          double dampingValue = integratedValue * velocity;

          damping->setValue(componentNo, dofNosLocal[i], dampingValue, ADD_VALUES);
        }
      }  // j
    }  // i
  }  // elementNo

  damping->finishGhostManipulation();

  LOG(DEBUG) << "damping: " << *damping;
}

template<typename Term>
void DynamicHyperelasticitySolver<Term>::
computeRungeKutta4()
{
  PetscErrorCode ierr;

  // compute k0 = Δt*a(u^(n), v^(n))
  computeAcceleration(u_, v_, k_[0]);
  ierr = VecScale(k_[0]->valuesGlobal(), this->timeStepWidth_); CHKERRV(ierr);

  // compute l0 = Δt*v^(n)
  ierr = VecCopy(v_->valuesGlobal(), l_[0]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecScale(l_[0]->valuesGlobal(), this->timeStepWidth_); CHKERRV(ierr);

  // compute k1 = Δt*a(u^(n) + 1/2*l0, v^(n) + 1/2*k0)
  // temp2 = u^(n) + 1/2*l0
  ierr = VecCopy(u_->valuesGlobal(), temp_[2]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(temp_[2]->valuesGlobal(), 0.5, l_[0]->valuesGlobal()); CHKERRV(ierr);

  // temp3 = v^(n) + 1/2*k0
  ierr = VecCopy(v_->valuesGlobal(), temp_[3]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(temp_[3]->valuesGlobal(), 0.5, k_[0]->valuesGlobal()); CHKERRV(ierr);
  computeAcceleration(temp_[2], temp_[3], k_[1]);

  // compute l1 = Δt*(v^(n) + 1/2*k0)
  ierr = VecCopy(temp_[3]->valuesGlobal(), l_[1]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecScale(l_[1]->valuesGlobal(), this->timeStepWidth_); CHKERRV(ierr);

  // compute k2 = Δt*a(u^(n) + 1/2*l1, v^(n) + 1/2*k1)
  // temp4 = u^(n) + 1/2*l1
  ierr = VecCopy(u_->valuesGlobal(), temp_[4]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(temp_[4]->valuesGlobal(), 0.5, l_[1]->valuesGlobal()); CHKERRV(ierr);

  // temp5 = v^(n) + 1/2*k1
  ierr = VecCopy(v_->valuesGlobal(), temp_[5]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(temp_[5]->valuesGlobal(), 0.5, k_[1]->valuesGlobal()); CHKERRV(ierr);
  computeAcceleration(temp_[4], temp_[5], k_[2]);

  // compute l2 = Δt*(v^(n) + 1/2*k1)
  ierr = VecCopy(temp_[5]->valuesGlobal(), l_[2]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecScale(l_[2]->valuesGlobal(), this->timeStepWidth_); CHKERRV(ierr);

  // compute k3 = Δt*a(u^(n) + l2, v^(n) + k2)
  // temp6 = u^(n) + l2
  ierr = VecCopy(u_->valuesGlobal(), temp_[6]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(temp_[6]->valuesGlobal(), 1.0, l_[2]->valuesGlobal()); CHKERRV(ierr);

  // temp7 = v^(n) + k2
  ierr = VecCopy(v_->valuesGlobal(), temp_[7]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(temp_[7]->valuesGlobal(), 1.0, k_[2]->valuesGlobal()); CHKERRV(ierr);
  computeAcceleration(temp_[6], temp_[7], k_[3]);

  // compute l3 = Δt*(v^(n) + k2)
  ierr = VecCopy(temp_[7]->valuesGlobal(), l_[3]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecScale(l_[3]->valuesGlobal(), this->timeStepWidth_); CHKERRV(ierr);

  // compute v^(n+1) = v^(n) + 1/6*(k0 + 2*k1 + 2*k2 + k3)
  ierr = VecAXPY(v_->valuesGlobal(), 1./6, k_[0]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(v_->valuesGlobal(), 2./6, k_[1]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(v_->valuesGlobal(), 2./6, k_[2]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(v_->valuesGlobal(), 1./6, k_[3]->valuesGlobal()); CHKERRV(ierr);

  // compute u^(n+1) = u^(n) + 1/6*(l0 + 2*l1 + 2*l2 + l3)
  ierr = VecAXPY(u_->valuesGlobal(), 1./6, l_[0]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(u_->valuesGlobal(), 2./6, l_[1]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(u_->valuesGlobal(), 2./6, l_[2]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecAXPY(u_->valuesGlobal(), 1./6, l_[3]->valuesGlobal()); CHKERRV(ierr);

  LOG(DEBUG) << "enforce incompressibility";

  // enforce incompressibility, solve K*u^(new) = K*u
  staticSolver_.materialComputeInternalVirtualWork(u_, temp_[8]);

  // solve ∂W_int(u) - temp8 = 0 with J = 1 for u. The previous values in u_ are the initial guess.
  staticSolver_.solveForDisplacements(temp_[8], u_);
}

template<typename Term>
void DynamicHyperelasticitySolver<Term>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

} // namespace TimeSteppingScheme

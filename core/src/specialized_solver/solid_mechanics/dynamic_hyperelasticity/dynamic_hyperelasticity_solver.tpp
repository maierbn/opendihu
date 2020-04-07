#include "specialized_solver/solid_mechanics/dynamic_hyperelasticity/dynamic_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.h"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

namespace TimeSteppingScheme
{

template<typename Term,typename MeshType>
DynamicHyperelasticitySolver<Term,MeshType>::
DynamicHyperelasticitySolver(DihuContext context) :
  TimeSteppingScheme(context["DynamicHyperelasticitySolver"]), hyperelasticitySolver_(context, "DynamicHyperelasticitySolver"), data_(context_)
{
  specificSettings_ = this->context_.getPythonConfig();
  density_ = specificSettings_.getOptionDouble("density", 1.0, PythonUtility::Positive);
  inputMeshIsGlobal_ = specificSettings_.getOptionBool("inputMeshIsGlobal", true);
  //viscosity_ = specificSettings_.getOptionDouble("viscosity", 0.0, PythonUtility::NonNegative);

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_["dynamic"], this->context_["dynamic"].getPythonConfig());
}

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
initialize()
{
  TimeSteppingScheme::initialize();

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("DynamicHyperelasticitySolver");

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  hyperelasticitySolver_.initialize();
  data_.setFunctionSpace(hyperelasticitySolver_.data().displacementsFunctionSpace());
  data_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // create variable of unknowns
  uvp_ = hyperelasticitySolver_.createPartitionedPetscVec("uvp");
  //uvp_ = hyperelasticitySolver_.combinedVecSolution();
  setInitialValues();

  PetscErrorCode ierr;
  ierr = VecDuplicate(uvp_->valuesGlobal(), &internalVirtualWork_); CHKERRV(ierr);
  ierr = VecDuplicate(uvp_->valuesGlobal(), &accelerationTerm_); CHKERRV(ierr);
  ierr = VecDuplicate(uvp_->valuesGlobal(), &externalVirtualWorkDead_); CHKERRV(ierr);

  LOG(DEBUG) << "internalVirtualWork_: " << internalVirtualWork_;
  LOG(DEBUG) << "accelerationTerm_: " << accelerationTerm_;
  LOG(DEBUG) << "externalVirtualWorkDead_: " << externalVirtualWorkDead_;
  LOG(DEBUG) << "uvp_: " << uvp_->valuesGlobal();

  // parse updateDirichletBoundaryConditionsFunction
  if (this->specificSettings_.hasKey("updateDirichletBoundaryConditionsFunction"))
  {
    PyObject *object = this->specificSettings_.getOptionPyObject("updateDirichletBoundaryConditionsFunction");
    if (object == Py_None)
    {
      pythonUpdateDirichletBoundaryConditionsFunction_ = nullptr;
    }
    else
    {
      pythonUpdateDirichletBoundaryConditionsFunction_ = this->specificSettings_.getOptionFunction("updateDirichletBoundaryConditionsFunction");
      updateDirichletBoundaryConditionsFunctionCallInterval_ = this->specificSettings_.getOptionInt("updateDirichletBoundaryConditionsFunctionCallInterval", 1, PythonUtility::Positive);
    }
  }

  // write initial mesh
  this->outputWriterManager_.writeOutput(this->data_, 0, 0.0);

  // check if initial values satisfy the static equation (for debugging)
  //hyperelasticitySolver_.debug();
}

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
callUpdateDirichletBoundaryConditionsFunction(double t)
{
  if (pythonUpdateDirichletBoundaryConditionsFunction_ == NULL)
    return;

  // only call this function at defined intervals
  if (updateDirichletBoundaryConditionsFunctionCallCount_ % updateDirichletBoundaryConditionsFunctionCallInterval_ != 0)
  {
    updateDirichletBoundaryConditionsFunctionCallCount_++;
    return;
  }
  updateDirichletBoundaryConditionsFunctionCallCount_++;

  // compose callback function
  PyObject *arglist = Py_BuildValue("(d)", t);
  PyObject *returnValue = PyObject_CallObject(pythonUpdateDirichletBoundaryConditionsFunction_, arglist);

  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();

  // parse the return value of the function
  std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBCValues =
    PythonUtility::convertFromPython<std::vector<std::pair<global_no_t,std::array<double,6>>>>::get(returnValue);

  LOG(DEBUG) << "newDirichletBCValues: " << newDirichletBCValues << "vecs: " << hyperelasticitySolver_.combinedVecSolution() << "," << uvp_;

  // set the new DirichletBC values
  hyperelasticitySolver_.combinedVecSolution()->updateDirichletBoundaryConditions(newDirichletBCValues, inputMeshIsGlobal_);
  uvp_->updateDirichletBoundaryConditions(newDirichletBCValues, inputMeshIsGlobal_);

  // decrement reference counters for python objects
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
setInitialValues()
{
  // set initial values as given in settings, or set to zero if not given
  std::vector<Vec3> localValuesDisplacements;
  std::vector<Vec3> localValuesVelocities;

  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace = this->data_.functionSpace();

  // determine if the initial values are given as global and local array
  bool inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {

    // if the settings specify a global list of values, extract the local values
    assert(displacementsFunctionSpace);

    // get number of global dofs, i.e. number of values in global list
    const int nDofsGlobal = displacementsFunctionSpace->nDofsGlobal();
    LOG(DEBUG) << "setInitialValues, nDofsGlobal = " << nDofsGlobal;

    // extract only the local dofs out of the list of global values
    this->specificSettings_.template getOptionVector<Vec3>("initialValuesDisplacements", nDofsGlobal, localValuesDisplacements);
    displacementsFunctionSpace->meshPartition()->extractLocalDofsWithoutGhosts(localValuesDisplacements);

    this->specificSettings_.template getOptionVector<Vec3>("initialValuesVelocities", nDofsGlobal, localValuesVelocities);
    displacementsFunctionSpace->meshPartition()->extractLocalDofsWithoutGhosts(localValuesVelocities);
  }
  else
  {
    // input is already only the local dofs, use all
    const int nDofsLocal = displacementsFunctionSpace->nDofsLocalWithoutGhosts();
    this->specificSettings_.template getOptionVector<Vec3>("initialValuesDisplacements", nDofsLocal, localValuesDisplacements);
    this->specificSettings_.template getOptionVector<Vec3>("initialValuesVelocities", nDofsLocal, localValuesVelocities);
  }
  VLOG(1) << "set initial values for displacements to " << localValuesDisplacements;
  VLOG(1) << "set initial values for velocities to " << localValuesVelocities;

  // set the first component of the solution variable by the given values
  this->data_.displacements()->setValuesWithoutGhosts(localValuesDisplacements);
  this->data_.velocities()->setValuesWithoutGhosts(localValuesVelocities);

  // set displacement entries in uvp
  int nDofsLocalWithoutGhosts = displacementsFunctionSpace->nDofsLocalWithoutGhosts();
  std::vector<double> localValues(nDofsLocalWithoutGhosts);

  uvp_->zeroGhostBuffer();

  // set displacement entries in uvp
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    for (int entryNo = 0; entryNo < nDofsLocalWithoutGhosts; entryNo++)
    {
      localValues[entryNo] = localValuesDisplacements[entryNo][componentNo];
    }

    uvp_->setValues(componentNo, nDofsLocalWithoutGhosts, displacementsFunctionSpace->meshPartition()->dofNosLocal().data(), localValues.data());
  }

  // set velocity entries in uvp
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    for (int entryNo = 0; entryNo < nDofsLocalWithoutGhosts; entryNo++)
    {
      localValues[entryNo] = localValuesVelocities[entryNo][componentNo];
    }

    uvp_->setValues(3+componentNo, nDofsLocalWithoutGhosts, displacementsFunctionSpace->meshPartition()->dofNosLocal().data(), localValues.data());
  }
  uvp_->zeroGhostBuffer();
  uvp_->finishGhostManipulation();
  //uvp_->setRepresentationGlobal();
}

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "DynamicHyperelasticitySolver::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "DynamicHyperelasticitySolver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    LOG(DEBUG) << "solveDynamicProblem";

    // advance simulation time, staticSolver needs the new time in order to output results after the solve is done
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // set the current Time to the hyperelasticity solver and then solve the dynamic problem
    hyperelasticitySolver_.setTimeSpan(-1, currentTime);
    hyperelasticitySolver_.solveDynamicProblem(uvp_, timeStepNo==0,
                                               internalVirtualWork_, externalVirtualWorkDead_, accelerationTerm_);

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // copy resulting values to data object such that they can be written by output writer
    hyperelasticitySolver_.setDisplacementsVelocitiesAndPressureFromCombinedVec(
      uvp_->valuesGlobal(), this->data_.displacements(), this->data_.velocities());

    hyperelasticitySolver_.setDisplacementsVelocitiesAndPressureFromCombinedVec(internalVirtualWork_, this->data_.internalVirtualWork());
    hyperelasticitySolver_.setDisplacementsVelocitiesAndPressureFromCombinedVec(externalVirtualWorkDead_, this->data_.externalVirtualWorkDead());
    hyperelasticitySolver_.setDisplacementsVelocitiesAndPressureFromCombinedVec(accelerationTerm_, this->data_.accelerationTerm());

    // write current output values
    this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);

    // potentially update DirichletBC by calling "updateDirichletBoundaryConditionsFunction"
    callUpdateDirichletBoundaryConditionsFunction(currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
    //this->data_.print();
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

template<typename Term,typename MeshType>
typename DynamicHyperelasticitySolver<Term,MeshType>::Data &DynamicHyperelasticitySolver<Term,MeshType>::
data()
{
  return data_;
}

template<typename Term,typename MeshType>
typename DynamicHyperelasticitySolver<Term,MeshType>::HyperelasticitySolverType &DynamicHyperelasticitySolver<Term,MeshType>::
hyperelasticitySolver()
{
  return hyperelasticitySolver_;
}


/*

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
initializeMassMatrix()
{
  LOG(TRACE) << "init mass matrix";

  // create inverseLumpedMassMatrix
  massMatrix_ = hyperelasticitySolver_.createPartitionedPetscMat("massMatrix");

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

  inverseLumpedMassMatrix_ = hyperelasticitySolver_.createPartitionedPetscMat("inverseLumpedMassMatrix");

  // row sum vector
  std::shared_ptr<VecHyperelasticity> rowSum = hyperelasticitySolver_.createPartitionedPetscVec("rowSum");

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

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
computeAcceleration(std::shared_ptr<VecHyperelasticity> u, std::shared_ptr<VecHyperelasticity> v, std::shared_ptr<VecHyperelasticity> a)
{
  assert(u);
  assert(v);
  assert(a);
  LOG(TRACE) << "computeAcceleration";

  // compute δW_int, store in temp_[0]
  hyperelasticitySolver_.materialComputeInternalVirtualWork(u, temp_[0]);

  LOG(DEBUG) << "δW_int: " << *temp_[0];
  // subtract external virtual work f

  // δW_ext = int_∂Ω T_a phi_L dS was precomputed in initialize of hyperelasticitySolver_
  PetscErrorCode ierr;
  ierr = VecAXPY(temp_[0]->valuesGlobal(), -1, hyperelasticitySolver_.externalVirtualWork()); CHKERRV(ierr);

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

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
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
        = SpatialDiscretization::IntegrandMassMatrix<D,EvaluationsType,DisplacementsFunctionSpace,1,double_v_t,dof_no_v_t,Term>::evaluateIntegrand(jacobian,xi) * viscosity_;

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

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
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

  // enforce incompressibility by solving K*u^(new) = K*u

  // compute temp8 = K*u
  //hyperelasticitySolver_.materialComputeInternalVirtualWork(u_, temp_[8]);

  // solve ∂W_int(u) - temp8 = 0 with J = 1 for u. The previous values in u_ are the initial guess.
  //hyperelasticitySolver_.solveForDisplacements(temp_[8], u_);
}

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
computeExplicitEuler()
{
  PetscErrorCode ierr;

  LOG(DEBUG) << "computeExplicitEuler, timeStepWidth_: " << this->timeStepWidth_;

  // compute k0 = Δt*a(u^(n), v^(n))
  computeAcceleration(u_, v_, k_[0]);
  LOG(DEBUG) << "k0: " << *k_[0];
  hyperelasticitySolver_.setDisplacementsAndPressureFromCombinedVec(k_[0]->valuesGlobal(), this->data_.acceleration());

  ierr = VecScale(k_[0]->valuesGlobal(), this->timeStepWidth_); CHKERRV(ierr);

  // compute l0 = Δt*v^(n)
  ierr = VecCopy(v_->valuesGlobal(), l_[0]->valuesGlobal()); CHKERRV(ierr);
  ierr = VecScale(l_[0]->valuesGlobal(), this->timeStepWidth_); CHKERRV(ierr);

  // compute v^(n+1) = v^(n) + k0
  ierr = VecAXPY(v_->valuesGlobal(), 1., k_[0]->valuesGlobal()); CHKERRV(ierr);

  // compute u^(n+1) = u^(n) + l0
  ierr = VecAXPY(u_->valuesGlobal(), 1., l_[0]->valuesGlobal()); CHKERRV(ierr);


  // copy the u^(n+1) and v^(n+1) values to data.rhs to be output after this time step
  hyperelasticitySolver_.setDisplacementsAndPressureFromCombinedVec(u_->valuesGlobal(), this->data_.displacementsCompressible());
  hyperelasticitySolver_.setDisplacementsAndPressureFromCombinedVec(v_->valuesGlobal(), this->data_.velocityCompressible());

  LOG(DEBUG) << "enforce incompressibility";

  // enforce incompressibility by solving K*u^(new) = K*u


  hyperelasticitySolver_.materialComputeInternalVirtualWork(u_, temp_[8]);
  hyperelasticitySolver_.solveForDisplacements(temp_[8], u_);    // solveForDisplacements(externalVirtualWork, displacements)
  return;
  // compute temp8 = K*u
  computeAcceleration(u_, v_, k_[0]);

  // temp0 = M*a^(n+1)
  ierr = MatMult(massMatrix_->valuesGlobal(), k_[0]->valuesGlobal(), temp_[0]->valuesGlobal()); CHKERRV(ierr);

  // δW_ext = int_∂Ω T_a phi_L dS was precomputed in initialize of hyperelasticitySolver_
  ierr = VecCopy(hyperelasticitySolver_.externalVirtualWork(), temp_[1]->valuesGlobal()); CHKERRV(ierr);

  // temp1 = temp1 - M*a^(n+1) = δW_ext - M*a^(n+1)
  VecAXPY(temp_[1]->valuesGlobal(), -1, temp_[0]->valuesGlobal());



  // copy the temp8 values to data.rhs to be output after this time step
  hyperelasticitySolver_.setDisplacementsAndPressureFromCombinedVec(temp_[1]->valuesGlobal(), this->data_.internalVirtualWorkCompressible());



  // solve ∂W_int(u) - temp8 = 0 with J = 1 for u. The previous values in u_ are the initial guess.
  hyperelasticitySolver_.solveForDisplacements(temp_[1], u_);    // solveForDisplacements(externalVirtualWork, displacements)

}*/


} // namespace TimeSteppingScheme

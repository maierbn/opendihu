#include "specialized_solver/solid_mechanics/dynamic_hyperelasticity/dynamic_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.h"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"
#include "spatial_discretization/neumann_boundary_conditions/01_neumann_boundary_conditions.h"

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
  pythonUpdateDirichletBoundaryConditionsFunction_ = nullptr;
  if (this->specificSettings_.hasKey("updateDirichletBoundaryConditionsFunction"))
  {
    PyObject *object = this->specificSettings_.getOptionPyObject("updateDirichletBoundaryConditionsFunction");
    if (object != Py_None)
    {
      pythonUpdateDirichletBoundaryConditionsFunction_ = this->specificSettings_.getOptionFunction("updateDirichletBoundaryConditionsFunction");
      updateDirichletBoundaryConditionsFunctionCallInterval_ = this->specificSettings_.getOptionInt("updateDirichletBoundaryConditionsFunctionCallInterval", 1, PythonUtility::Positive);
    }
  }

  // parse updateNeumannBoundaryConditionsFunction
  pythonUpdateNeumannBoundaryConditionsFunction_ = nullptr;
  if (this->specificSettings_.hasKey("updateNeumannBoundaryConditionsFunction"))
  {
    LOG(DEBUG) << "parse updateNeumannBoundaryConditionsFunction";
    PyObject *object = this->specificSettings_.getOptionPyObject("updateNeumannBoundaryConditionsFunction");
    if (object != Py_None)
    {
      pythonUpdateNeumannBoundaryConditionsFunction_ = this->specificSettings_.getOptionFunction("updateNeumannBoundaryConditionsFunction");
      updateNeumannBoundaryConditionsFunctionCallInterval_ = this->specificSettings_.getOptionInt("updateNeumannBoundaryConditionsFunctionCallInterval", 1, PythonUtility::Positive);
      LOG(DEBUG) << "interval: " << updateNeumannBoundaryConditionsFunctionCallInterval_;
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

  updateDirichletBoundaryConditions(newDirichletBCValues);

  // decrement reference counters for python objects
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
updateDirichletBoundaryConditions(std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBCValues)
{
  // set the new DirichletBC values
  hyperelasticitySolver_.combinedVecSolution()->updateDirichletBoundaryConditions(newDirichletBCValues, inputMeshIsGlobal_);
  uvp_->updateDirichletBoundaryConditions(newDirichletBCValues, inputMeshIsGlobal_);
}

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
addDirichletBoundaryConditions(std::vector<typename SpatialDiscretization::DirichletBoundaryConditions<DisplacementsFunctionSpace,6>::ElementWithNodes> &boundaryConditionElements, bool overwriteBcOnSameDof)
{
  hyperelasticitySolver_.addDirichletBoundaryConditions(boundaryConditionElements, overwriteBcOnSameDof);

  // recreate all vectors
  uvp_ = hyperelasticitySolver_.createPartitionedPetscVec("uvp");

  //uvp_->startGhostManipulation();
  uvp_->zeroGhostBuffer();
  uvp_->finishGhostManipulation();

  PetscErrorCode ierr;
  ierr = VecDuplicate(uvp_->valuesGlobal(), &internalVirtualWork_); CHKERRV(ierr);
  ierr = VecDuplicate(uvp_->valuesGlobal(), &accelerationTerm_); CHKERRV(ierr);
  ierr = VecDuplicate(uvp_->valuesGlobal(), &externalVirtualWorkDead_); CHKERRV(ierr);

}

template<typename Term,typename MeshType>
void DynamicHyperelasticitySolver<Term,MeshType>::
callUpdateNeumannBoundaryConditionsFunction(double t)
{
  if (pythonUpdateNeumannBoundaryConditionsFunction_ == NULL)
    return;

  // only call this function at defined intervals
  if (updateNeumannBoundaryConditionsFunctionCallCount_ % updateNeumannBoundaryConditionsFunctionCallInterval_ != 0)
  {
    updateNeumannBoundaryConditionsFunctionCallCount_++;
    return;
  }
  updateNeumannBoundaryConditionsFunctionCallCount_++;

  // compose callback function
  PyObject *arglist = Py_BuildValue("(d)", t);
  PyObject *returnValue = PyObject_CallObject(pythonUpdateNeumannBoundaryConditionsFunction_, arglist);

  PythonUtility::checkForError();

  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();

  // parse the return value of the function
  using NeumannBoundaryConditionsType = typename SpatialDiscretization::NeumannBoundaryConditions<typename HyperelasticitySolverType::DisplacementsFunctionSpace,Quadrature::Gauss<3>,3>;
  std::shared_ptr<NeumannBoundaryConditionsType> newNeumannBoundaryConditions = std::make_shared<NeumannBoundaryConditionsType>(this->context_);
  newNeumannBoundaryConditions->initialize(returnValue, this->data_.functionSpace(), "neumannBoundaryConditions");

  hyperelasticitySolver_.updateNeumannBoundaryConditions(newNeumannBoundaryConditions);

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
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
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

    // potentially update NeumannBC by calling "updateNeumannBoundaryConditionsFunction"
    callUpdateNeumannBoundaryConditionsFunction(currentTime);

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

} // namespace TimeSteppingScheme

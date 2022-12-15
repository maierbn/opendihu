#include "specialized_solver/solid_mechanics/quasistatic_hyperelasticity/quasistatic_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.h"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"
#include "spatial_discretization/neumann_boundary_conditions/01_neumann_boundary_conditions.h"

namespace TimeSteppingScheme
{

template<typename Term,bool withLargeOutput,typename MeshType>
QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
QuasistaticHyperelasticitySolver(DihuContext context) :
  TimeSteppingScheme(context["QuasistaticHyperelasticitySolver"]), hyperelasticitySolver_(context, "QuasistaticHyperelasticitySolver"), data_(context_)
{
  specificSettings_ = this->context_.getPythonConfig();
  density_ = specificSettings_.getOptionDouble("density", 1.0, PythonUtility::Positive);
  inputMeshIsGlobal_ = specificSettings_.getOptionBool("inputMeshIsGlobal", true);
  //viscosity_ = specificSettings_.getOptionDouble("viscosity", 0.0, PythonUtility::NonNegative);

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_["dynamic"], this->context_["dynamic"].getPythonConfig());
  isReferenceGeometryInitialized_ = false;
}

template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
initialize()
{
  LOG_SCOPE_FUNCTION;
  TimeSteppingScheme::initialize();

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("QuasistaticHyperelasticitySolver");

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  hyperelasticitySolver_.initialize();
  data_.setFunctionSpace(hyperelasticitySolver_.data().displacementsFunctionSpace());
  data_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // determine if there is any traction defined in current configuration, if so it has to be converted to reference configuration in every timestep
  isTractionInCurrentConfiguration_ = hyperelasticitySolver_.neumannBoundaryConditions()->isTractionInCurrentConfiguration();

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

  pythonTotalForceFunction_ = nullptr;
  if (this->specificSettings_.hasKey("totalForceFunction"))
  {
    LOG(DEBUG) << "parse totalForceFunction";
    PyObject *object = this->specificSettings_.getOptionPyObject("totalForceFunction");
    if (object != Py_None)
    {
      pythonTotalForceFunction_ = this->specificSettings_.getOptionFunction("totalForceFunction");
      pythonTotalForceFunctionCallInterval_ = this->specificSettings_.getOptionInt("totalForceFunctionCallInterval", 1, PythonUtility::Positive);
    }
  }

  // parse elements that are used to compute total force at top and bottom
  if (this->specificSettings_.hasKey("totalForceTopElementNosGlobal"))
  {
    std::vector<global_no_t> totalForceTopElementsGlobal;
    this->specificSettings_.getOptionVector<global_no_t>("totalForceTopElementNosGlobal", totalForceTopElementsGlobal);
    std::set<global_no_t> totalForceTopElementsGlobalSet(totalForceTopElementsGlobal.begin(), totalForceTopElementsGlobal.end());

    LOG(DEBUG) << "given elements top: " << totalForceTopElementsGlobalSet;

    // iterate over local elements
    for (element_no_t elementNoLocal = 0; elementNoLocal < this->data_.functionSpace()->nElementsLocal(); elementNoLocal++)
    {
      global_no_t elementNoGlobal = this->data_.functionSpace()->meshPartition()->getElementNoGlobalNatural(elementNoLocal);

      if (totalForceTopElementsGlobalSet.find(elementNoGlobal) != totalForceTopElementsGlobalSet.end())
      {
        LOG(DEBUG) << " (" << elementNoLocal << "->" << elementNoGlobal << ") yes";
        bottomTopElements_.push_back(std::tuple<element_no_t,bool>(elementNoLocal,true));
      }
      else
        LOG(DEBUG) << " (" << elementNoLocal << "->" << elementNoGlobal << ") no";
    }
  }

  if (this->specificSettings_.hasKey("totalForceBottomElementNosGlobal"))
  {
    std::vector<global_no_t> totalForceBottomElementsGlobal;
    this->specificSettings_.getOptionVector<global_no_t>("totalForceBottomElementNosGlobal", totalForceBottomElementsGlobal);
    std::set<global_no_t> totalForceBottomElementsGlobalSet(totalForceBottomElementsGlobal.begin(), totalForceBottomElementsGlobal.end());

    LOG(DEBUG) << "given elements bottom: " << totalForceBottomElementsGlobalSet;
    // iterate over local elements
    for (element_no_t elementNoLocal = 0; elementNoLocal < this->data_.functionSpace()->nElementsLocal(); elementNoLocal++)
    {
      global_no_t elementNoGlobal = this->data_.functionSpace()->meshPartition()->getElementNoGlobalNatural(elementNoLocal);

      if (totalForceBottomElementsGlobalSet.find(elementNoGlobal) != totalForceBottomElementsGlobalSet.end())
      {
        LOG(DEBUG) << " (" << elementNoLocal << "->" << elementNoGlobal << ") yes";
        bottomTopElements_.push_back(std::tuple<element_no_t,bool>(elementNoLocal,false));
      }
      else
        LOG(DEBUG) << " (" << elementNoLocal << "->" << elementNoGlobal << ") no";
    }
  }

  totalForceLogFilename_ = this->specificSettings_.getOptionString("totalForceLogFilename", "");

  // write initial mesh but don't increment counter
  this->outputWriterManager_.writeOutput(this->data_, 0, 0.0, 0);

  // check if initial values satisfy the static equation (for debugging)
  //hyperelasticitySolver_.debug();
}

template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
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

template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
updateDirichletBoundaryConditions(std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBCValues)
{
  // set the new DirichletBC values
  hyperelasticitySolver_.combinedVecSolution()->updateDirichletBoundaryConditions(newDirichletBCValues, inputMeshIsGlobal_);
  uvp_->updateDirichletBoundaryConditions(newDirichletBCValues, inputMeshIsGlobal_);
}

template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
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

template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
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
  newNeumannBoundaryConditions->setDeformationGradientField(this->hyperelasticitySolver_.data().deformationGradient());

  isTractionInCurrentConfiguration_ = newNeumannBoundaryConditions->isTractionInCurrentConfiguration();

  hyperelasticitySolver_.updateNeumannBoundaryConditions(newNeumannBoundaryConditions);

  // decrement reference counters for python objects
  Py_CLEAR(returnValue);
  Py_CLEAR(arglist);
}

template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
updateNeumannBoundaryConditions()
{
  LOG(INFO) << "Updating right hand side because traction was specified in current configuration.";

  // update rhs with new traction values in reference configuration if there are moving loads (traction in current configuration)
  hyperelasticitySolver_.neumannBoundaryConditions()->initializeRhs();

  // update δW_ext,dead
  hyperelasticitySolver_.updateNeumannBoundaryConditions(hyperelasticitySolver_.neumannBoundaryConditions());
}

template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
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

template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
advanceTimeSpan(bool withOutputWritersEnabled)
{
  LOG_SCOPE_FUNCTION;
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // update the reference geometry to use the current value, in the first time step
  if (!isReferenceGeometryInitialized_)
  {
    this->hyperelasticitySolver_.data().updateReferenceGeometry();
    isReferenceGeometryInitialized_ = true;
  }

  // only before the first timestep
  // in case there have been new displacement values set by the data connector slots, update the uvp_ variable
  if (this->startTime_ < 1e-12)
  {
    // loop over x,y,z components of displacements
    static std::vector<double> values;
    for (int i = 0; i < 3; i++)
    {
      values.clear();
      this->data_.displacements()->getValuesWithoutGhosts(i, values);

      LOG(DEBUG) << "component " << i << ", initialize displacement values: " << values;
      uvp_->setValues(i, data_.functionSpace()->meshPartition()->nDofsLocalWithoutGhosts(),
                      data_.functionSpace()->meshPartition()->dofNosLocal().data(), values.data(), INSERT_VALUES);
    }
    // representation is combined-global
    // assemble vector
    PetscErrorCode ierr;
    ierr = VecAssemblyBegin(uvp_->valuesGlobal()); CHKERRV(ierr);
    ierr = VecAssemblyEnd(uvp_->valuesGlobal()); CHKERRV(ierr);
  }

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "QuasistaticHyperelasticitySolver::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "QuasistaticHyperelasticitySolver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    LOG(DEBUG) << "solveStaticProblem";

    // advance simulation time, staticSolver needs the new time in order to output results after the solve is done
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // set the current Time to the hyperelasticity solver and then solve the dynamic problem
    hyperelasticitySolver_.setTimeSpan(-1, currentTime);
    hyperelasticitySolver_.solveQuasistaticProblem(uvp_, timeStepNo==0, withOutputWritersEnabled);    // stop duration measurement
    
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // copy resulting values to data object such that they can be written by output writer
    hyperelasticitySolver_.setDisplacementsVelocitiesAndPressureFromCombinedVec(
     uvp_->valuesGlobal(), this->data_.displacements(), this->data_.velocities());

    hyperelasticitySolver_.setDisplacementsVelocitiesAndPressureFromCombinedVec(internalVirtualWork_, this->data_.internalVirtualWork());
    hyperelasticitySolver_.setDisplacementsVelocitiesAndPressureFromCombinedVec(externalVirtualWorkDead_, this->data_.externalVirtualWorkDead());
    hyperelasticitySolver_.setDisplacementsVelocitiesAndPressureFromCombinedVec(accelerationTerm_, this->data_.accelerationTerm());

  
    if (withOutputWritersEnabled)
      this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);

      // this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);
        
    // potentially update DirichletBC by calling "updateDirichletBoundaryConditionsFunction"
    callUpdateDirichletBoundaryConditionsFunction(currentTime);

    // potentially update NeumannBC by calling "updateNeumannBoundaryConditionsFunction"
    callUpdateNeumannBoundaryConditionsFunction(currentTime);

    if (isTractionInCurrentConfiguration_)
      updateNeumannBoundaryConditions();

    // compute the total force and torque at the z+ and z- surfaces of the volume
    computeBearingForcesAndMoments(currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
    //this->data_.print();
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

//! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  //call the output writer of the nested solver
  this->hyperelasticitySolver_.callOutputWriter(timeStepNo, currentTime, callCountIncrement);

  // call the own output writer
  this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime, callCountIncrement);
}

template<typename Term,bool withLargeOutput,typename MeshType>
void QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
computeBearingForcesAndMoments(double currentTime)
{
  if (totalForceLogFilename_.empty())
    return;

  LOG(DEBUG) << bottomTopElements_.size() << " elements: " << bottomTopElements_;

  // compute the total forces and moments
  Vec3 bearingForceBottom;
  Vec3 bearingMomentBottom;
  Vec3 bearingForceTop;
  Vec3 bearingMomentTop;
  hyperelasticitySolver_.computeBearingForceAndMoment(bottomTopElements_,
    bearingForceBottom, bearingMomentBottom, bearingForceTop, bearingMomentTop);

  if (DihuContext::ownRankNoCommWorld() == 0)
  {
    static bool headerWritten = false;

    // append to file
    std::ofstream file;
    OutputWriter::Generic::openFile(file, totalForceLogFilename_, true);

    if (!headerWritten)
    {
      file << "currentTime;forceBottomX;forceBottomY;forceBottomZ;momentBottomX;momentBottomY;momentBottomZ;"
        << "forceTopX;forceTopY;forceTopZ;momentTopX;momentTopY;momentTopZ\n";
      headerWritten = true;
    }

    // write currentTime;forceBottom;momentBottom;forceTop;momentTop
    file << currentTime;

    for (int i = 0; i < 3; i++)
      file << ";" << bearingForceBottom[i];

    for (int i = 0; i < 3; i++)
      file << ";" << bearingMomentBottom[i];

    for (int i = 0; i < 3; i++)
      file << ";" << bearingForceTop[i];

    for (int i = 0; i < 3; i++)
      file << ";" << bearingMomentTop[i];
    file << std::endl;
  }

  // call python callback
  if (pythonTotalForceFunction_)
  {
    // only call this function at defined intervals
    if (pythonTotalForceFunctionCallCount_ % pythonTotalForceFunctionCallInterval_ != 0)
    {
      pythonTotalForceFunctionCallCount_++;
      return;
    }
    pythonTotalForceFunctionCallCount_++;

    // create four python variables
    PyObject *bearingForceBottomList = PythonUtility::convertToPython<Vec3>::get(bearingForceBottom);
    PyObject *bearingMomentBottomList = PythonUtility::convertToPython<Vec3>::get(bearingMomentBottom);
    PyObject *bearingForceTopList = PythonUtility::convertToPython<Vec3>::get(bearingForceTop);
    PyObject *bearingMomentTopList = PythonUtility::convertToPython<Vec3>::get(bearingMomentTop);

    // compose callback function
    PyObject *arglist = Py_BuildValue("(d,O,O,O,O)", currentTime, bearingForceBottomList, bearingMomentBottomList,
                                      bearingForceTopList, bearingMomentTopList);
    PyObject *returnValue = PyObject_CallObject(pythonTotalForceFunction_, arglist);

    PythonUtility::checkForError();

    // if there was an error while executing the function, print the error message
    if (returnValue == NULL)
      PyErr_Print();
  }
}

template<typename Term,bool withLargeOutput,typename MeshType>
typename QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::Data &QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
data()
{
  return data_;
}

template<typename Term,bool withLargeOutput,typename MeshType>
Vec QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
currentState()
{
  return uvp_->valuesGlobal();
}

template<typename Term,bool withLargeOutput,typename MeshType>
typename QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::HyperelasticitySolverType &QuasistaticHyperelasticitySolver<Term,withLargeOutput,MeshType>::
hyperelasticitySolver()
{
  return hyperelasticitySolver_;
}

} // namespace TimeSteppingScheme

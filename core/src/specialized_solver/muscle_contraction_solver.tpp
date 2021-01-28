#include "specialized_solver/muscle_contraction_solver.h"

#include <omp.h>
#include <sstream>

#include "utility/math_utility.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

template<typename MeshType,typename Term,bool withLargeOutputFiles>
MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
MuscleContractionSolver(DihuContext context) :
  Runnable(),
  ::TimeSteppingScheme::TimeSteppingScheme(context["MuscleContractionSolver"]),
  data_(this->context_), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();

  // determine if the dynamic or the quasi-static formulation is used
  isDynamic_ = this->specificSettings_.getOptionBool("dynamic", true);

  if (isDynamic_)
  {
    LOG(DEBUG) << "use dynamic hyperelasticity solver";
    dynamicHyperelasticitySolver_ = std::make_shared<DynamicHyperelasticitySolverType>(this->context_);
  }
  else
  {
    LOG(DEBUG) << "use static hyperelasticity solver";
    staticHyperelasticitySolver_ = std::make_shared<StaticHyperelasticitySolverType>(this->context_);
  }

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  // parse options
  pmax_ = this->specificSettings_.getOptionDouble("Pmax", 1.0, PythonUtility::Positive);
  enableForceLengthRelation_ = this->specificSettings_.getOptionBool("enableForceLengthRelation", true);
  lambdaDotScalingFactor_ = this->specificSettings_.getOptionDouble("lambdaDotScalingFactor", 1.0);

  this->specificSettings_.template getOptionVector<std::string>("mapGeometryToMeshes", meshNamesOfGeometryToMapTo_);
}

template<typename MeshType,typename Term,bool withLargeOutputFiles>
void MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
advanceTimeSpan(bool withOutputWritersEnabled)
{
  // This method computes some time steps of the simulation by running a for loop over the time steps.
  // The number of steps, timestep width and current time are all set by the parent class, TimeSteppingScheme.

  // create mapping between meshes for the geometry field mapping, do this here already at the beginning where the meshes are not yet deformed
  initializeMappingBetweenMeshes();

  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute time span of this method
  double timeSpan = this->endTime_ - this->startTime_;

  LOG_N_TIMES(3,INFO) << "durationComputeMappingBetweenMeshes: " << Control::PerformanceMeasurement::getDuration("durationComputeMappingBetweenMeshes");

  // output for debugging
  LOG(DEBUG) << "MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    // in defined intervals (settings "timeStepOutputInterval") print out the current timestep
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "MuscleContractionSolver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // compute the current active stress
    computeActiveStress();

    if (isDynamic_)
    {
      this->dynamicHyperelasticitySolver_->setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

      LOG(DEBUG) << "call dynamic hyperelasticitySolver";

      // advance the simulation by the specified time span
      dynamicHyperelasticitySolver_->advanceTimeSpan(withOutputWritersEnabled);
    }
    else
    {
      staticHyperelasticitySolver_->run();
    }

    // compute new values of λ and λ_dot, to be transferred to the CellML solvers
    computeLambda();

    // advance simulation time
    timeStepNo++;

    // compute new current simulation time
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values using the output writers
    if (withOutputWritersEnabled)
      this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  mapGeometryToGivenMeshes();
}

template<typename MeshType,typename Term,bool withLargeOutputFiles>
void MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
initialize()
{
  // only initialize once
  if (initialized_)
    return;

  LOG(DEBUG) << "MuscleContractionSolver::initialize";

  // initialize() will be called before the simulation starts.

  // call initialize of the parent class, this parses the timestepping settings from the settings file
  TimeSteppingScheme::TimeSteppingScheme::initialize();

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("MuscleContractionSolver");

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested timestepping solver
  if (isDynamic_)
  {
    dynamicHyperelasticitySolver_->initialize();
  }
  else
  {
    staticHyperelasticitySolver_->initialize();
  }

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // In order to initialize the data object and actuall create all variables, we first need to assign a function space to the data object.
  // A function space object of type FunctionSpace<MeshType,BasisFunctionType> (see "function_space/function_space.h")
  // is an object that stores the mesh (e.g., all nodes and elements) as well as the basis function (e.g. linear Lagrange basis functions).
  // The dynamicHyperelasticitySolver_ solver already created a function space that we should use. We already have a typedef "FunctionSpace" that is the class of dynamicHyperelasticitySolver_'s function space type.
  std::shared_ptr<FunctionSpace> functionSpace;
  if (isDynamic_)
    functionSpace = dynamicHyperelasticitySolver_->data().functionSpace();
  else
    functionSpace = staticHyperelasticitySolver_->data().functionSpace();

  // Pass the function space to the data object. data_ stores field variables.
  // It needs to know the number of nodes and degrees of freedom (dof) from the function space in order to create the vectors with the right size.
  data_.setFunctionSpace(functionSpace);

  // now call initialize, data will then create all variables (Petsc Vec's)
  data_.initialize();

  bool setGeometryFieldForTransfer = meshNamesOfGeometryToMapTo_.empty();

  if (isDynamic_)
  {
    typename DynamicHyperelasticitySolverType::HyperelasticitySolverType &hyperelasticitySolver = dynamicHyperelasticitySolver_->hyperelasticitySolver();

    // set field variables from dynamicHyperelasticitySolver in data_ such that they can be output by the output writer
    data_.setFieldVariables(dynamicHyperelasticitySolver_->data().displacements(),
                            dynamicHyperelasticitySolver_->data().velocities(),
                            hyperelasticitySolver.data().activePK2Stress(),
                            hyperelasticitySolver.data().pK2Stress(),
                            hyperelasticitySolver.data().fiberDirection(),
                            hyperelasticitySolver.data().materialTraction(),
                            setGeometryFieldForTransfer);
  }
  else
  {
    // set field variables from dynamicHyperelasticitySolver in data_ such that they can be output by the output writer
    data_.setFieldVariables(staticHyperelasticitySolver_->data().displacements(),
                            staticHyperelasticitySolver_->data().velocities(),
                            staticHyperelasticitySolver_->data().activePK2Stress(),
                            staticHyperelasticitySolver_->data().pK2Stress(),
                            staticHyperelasticitySolver_->data().fiberDirection(),
                            staticHyperelasticitySolver_->data().materialTraction(),
                            setGeometryFieldForTransfer);
  }

  // set the slotConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setSlotConnectorData(getSlotConnectorData());

  initialized_ = true;
}

template<typename MeshType,typename Term,bool withLargeOutputFiles>
void MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
run()
{
  initialize();

  advanceTimeSpan();
}

template<typename MeshType,typename Term,bool withLargeOutputFiles>
void MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
reset()
{
  if (isDynamic_)
    dynamicHyperelasticitySolver_->reset();
  else
    staticHyperelasticitySolver_->reset();

  // "uninitialize" everything
}

//! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
template<typename MeshType,typename Term,bool withLargeOutputFiles>
void MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  // write output of nested solver
  if (isDynamic_)
  {
    this->dynamicHyperelasticitySolver_->callOutputWriter(timeStepNo, currentTime, callCountIncrement);
  }
  else
  {
    staticHyperelasticitySolver_->callOutputWriter(timeStepNo, currentTime, callCountIncrement);
  }

  // write output of own output writers
  this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime, callCountIncrement);
}

template<typename MeshType,typename Term,bool withLargeOutputFiles>
void MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
initializeMappingBetweenMeshes()
{
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_+std::string("_map_geometry"));

  LOG(INFO) << "initializeMappingBetweenMeshes, meshNamesOfGeometryToMapTo_=" << meshNamesOfGeometryToMapTo_;

  if (!meshNamesOfGeometryToMapTo_.empty())
  {
    using SourceFunctionSpaceType = typename StaticHyperelasticitySolverType::DisplacementsFunctionSpace;
    //using SourceFieldVariableType = FieldVariable::FieldVariable<SourceFunctionSpaceType,3>;

    assert(data_.functionSpace());

    // get source function space
    std::shared_ptr<SourceFunctionSpaceType> functionSpaceSource = data_.functionSpace();

    // loop over all given mesh names to which we should transfer the geometry
    for (std::string meshName : meshNamesOfGeometryToMapTo_)
    {

      // for first order meshes
      using TargetFunctionSpaceType1 = ::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>;
      LOG(DEBUG) << "mesh \"" << meshName << "\", test if " << StringUtility::demangle(typeid(TargetFunctionSpaceType1).name());

      // if the mesh name corresponds to a linear mesh
      if (this->context_.meshManager()->hasFunctionSpaceOfType<TargetFunctionSpaceType1>(meshName))
      {
        // get target function space
        std::shared_ptr<TargetFunctionSpaceType1> functionSpaceTarget = this->context_.meshManager()->functionSpace<TargetFunctionSpaceType1>(meshName);

        LOG(DEBUG) << "** create mapping " << functionSpaceSource->meshName() << " -> " << functionSpaceTarget->meshName();

        // create mapping between functionSpaceSource and functionSpaceTarget
        DihuContext::mappingBetweenMeshesManager()->template mappingBetweenMeshes<SourceFunctionSpaceType,TargetFunctionSpaceType1>(functionSpaceSource, functionSpaceTarget);
      }
      else LOG(DEBUG) << "no";

      // for second order meshes
      using TargetFunctionSpaceType2 = ::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>;
      LOG(DEBUG) << "mesh \"" << meshName << "\", test if " << StringUtility::demangle(typeid(TargetFunctionSpaceType2).name());

      // if the mesh name corresponds to a quadratic mesh
      if (this->context_.meshManager()->hasFunctionSpaceOfType<TargetFunctionSpaceType2>(meshName))
      {
        // get target function space
        std::shared_ptr<TargetFunctionSpaceType2> functionSpaceTarget = this->context_.meshManager()->functionSpace<TargetFunctionSpaceType2>(meshName);

        LOG(DEBUG) << "** create mapping " << functionSpaceSource->meshName() << " -> " << functionSpaceTarget->meshName();

        // create mapping between functionSpaceSource and functionSpaceTarget
        DihuContext::mappingBetweenMeshesManager()->template mappingBetweenMeshes<SourceFunctionSpaceType,TargetFunctionSpaceType2>(functionSpaceSource, functionSpaceTarget);
      }
      else LOG(DEBUG) << "no";

      // for first order composite meshes
      using TargetFunctionSpaceType3 = ::FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>;
      LOG(DEBUG) << "mesh \"" << meshName << "\", test if " << StringUtility::demangle(typeid(TargetFunctionSpaceType3).name());

      // if the mesh name corresponds to a linear mesh
      if (this->context_.meshManager()->hasFunctionSpaceOfType<TargetFunctionSpaceType3>(meshName))
      {
        // get target function space
        std::shared_ptr<TargetFunctionSpaceType3> functionSpaceTarget = this->context_.meshManager()->functionSpace<TargetFunctionSpaceType3>(meshName);

        LOG(DEBUG) << "** create mapping " << functionSpaceSource->meshName() << " -> " << functionSpaceTarget->meshName();

        // create mapping between functionSpaceSource and functionSpaceTarget
        //DihuContext::mappingBetweenMeshesManager()->template mappingBetweenMeshes<SourceFunctionSpaceType,TargetFunctionSpaceType3>(functionSpaceSource, functionSpaceTarget);
        DihuContext::mappingBetweenMeshesManager()->template mappingBetweenMeshes<TargetFunctionSpaceType3,SourceFunctionSpaceType>(functionSpaceTarget, functionSpaceSource);
      }
      else LOG(DEBUG) << "no";

      // for second order composite meshes
      using TargetFunctionSpaceType4 = ::FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>;
      LOG(DEBUG) << "mesh \"" << meshName << "\", test if " << StringUtility::demangle(typeid(TargetFunctionSpaceType4).name());

      // if the mesh name corresponds to a quadratic mesh
      if (this->context_.meshManager()->hasFunctionSpaceOfType<TargetFunctionSpaceType4>(meshName))
      {
        // get target function space
        std::shared_ptr<TargetFunctionSpaceType4> functionSpaceTarget = this->context_.meshManager()->functionSpace<TargetFunctionSpaceType4>(meshName);

        LOG(DEBUG) << "** create mapping " << functionSpaceSource->meshName() << " -> " << functionSpaceTarget->meshName();

        // create mapping between functionSpaceSource and functionSpaceTarget
        //DihuContext::mappingBetweenMeshesManager()->template mappingBetweenMeshes<SourceFunctionSpaceType,TargetFunctionSpaceType4>(functionSpaceSource, functionSpaceTarget);
        DihuContext::mappingBetweenMeshesManager()->template mappingBetweenMeshes<TargetFunctionSpaceType4,SourceFunctionSpaceType>(functionSpaceTarget, functionSpaceSource);
      }
      else LOG(DEBUG) << "no";
    }
  }

  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_+std::string("_map_geometry"));
}

template<typename MeshType,typename Term,bool withLargeOutputFiles>
void MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
mapGeometryToGivenMeshes()
{
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_+std::string("_map_geometry"));

  LOG(DEBUG) << "mapGeometryToGivenMeshes: meshNamesOfGeometryToMapTo: " << meshNamesOfGeometryToMapTo_;
  if (!meshNamesOfGeometryToMapTo_.empty())
  {
    using SourceFunctionSpaceType = typename StaticHyperelasticitySolverType::DisplacementsFunctionSpace;
    using SourceFieldVariableType = FieldVariable::FieldVariable<SourceFunctionSpaceType,3>;
    
    assert(data_.functionSpace());

    // get source field variable
    std::shared_ptr<SourceFieldVariableType> geometryFieldSource = std::make_shared<SourceFieldVariableType>(data_.functionSpace()->geometryField());

    std::vector<Vec3> geometryValuesSource;
    geometryFieldSource->getValuesWithoutGhosts(geometryValuesSource);

    LOG(DEBUG) << "geometryValuesSource: " << geometryValuesSource;

    // loop over all given mesh names to which we should transfer the geometry
    for (std::string meshName : meshNamesOfGeometryToMapTo_)
    {
      // for first order meshes
      using TargetFunctionSpaceType1 = ::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>;
      using TargetFieldVariableType1 = FieldVariable::FieldVariable<TargetFunctionSpaceType1,3>;

      // if the mesh name corresponds to a linear mesh
      if (this->context_.meshManager()->hasFunctionSpaceOfType<TargetFunctionSpaceType1>(meshName))
      {
        // get target geometry field variable
        std::shared_ptr<TargetFieldVariableType1> geometryFieldTarget = std::make_shared<TargetFieldVariableType1>(
          this->context_.meshManager()->functionSpace<TargetFunctionSpaceType1>(meshName)->geometryField());

        LOG(DEBUG) << "transfer geometry field to linear mesh, " << geometryFieldSource->functionSpace()->meshName() << " -> "
          << geometryFieldTarget->functionSpace()->meshName();
        LOG(DEBUG) << StringUtility::demangle(typeid(SourceFunctionSpaceType).name()) << " -> " << StringUtility::demangle(typeid(TargetFunctionSpaceType1).name());

        // perform the mapping
        DihuContext::mappingBetweenMeshesManager()->template prepareMapping<SourceFieldVariableType,TargetFieldVariableType1>(geometryFieldSource, geometryFieldTarget, -1);

        // map the whole geometry field (all components), do not avoid copy
        DihuContext::mappingBetweenMeshesManager()->template map<SourceFieldVariableType,TargetFieldVariableType1>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
        DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<SourceFieldVariableType,TargetFieldVariableType1>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
      }

      // for second order meshes
      using TargetFunctionSpaceType2 = ::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>;
      using TargetFieldVariableType2 = FieldVariable::FieldVariable<TargetFunctionSpaceType2,3>;

      // if the mesh name corresponds to a quadratic mesh
      if (this->context_.meshManager()->hasFunctionSpaceOfType<TargetFunctionSpaceType2>(meshName))
      {
        // get target geometry field variable
        std::shared_ptr<TargetFieldVariableType2> geometryFieldTarget = std::make_shared<TargetFieldVariableType2>(
          this->context_.meshManager()->functionSpace<TargetFunctionSpaceType2>(meshName)->geometryField());

        LOG(DEBUG) << "transfer geometry field to quadratic mesh, " << geometryFieldSource->functionSpace()->meshName() << " -> "
          << geometryFieldTarget->functionSpace()->meshName();
        LOG(DEBUG) << StringUtility::demangle(typeid(SourceFunctionSpaceType).name()) << " -> " << StringUtility::demangle(typeid(TargetFunctionSpaceType2).name());

        // perform the mapping
        DihuContext::mappingBetweenMeshesManager()->template prepareMapping<SourceFieldVariableType,TargetFieldVariableType2>(geometryFieldSource, geometryFieldTarget, -1);

        // map the whole geometry field (all components), do not avoid copy
        DihuContext::mappingBetweenMeshesManager()->template map<SourceFieldVariableType,TargetFieldVariableType2>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
        DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<SourceFieldVariableType,TargetFieldVariableType2>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
      }

      // for first order composite meshes
      using TargetFunctionSpaceType3 = ::FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>;
      using TargetFieldVariableType3 = FieldVariable::FieldVariable<TargetFunctionSpaceType3,3>;

      // if the mesh name corresponds to a linear mesh
      if (this->context_.meshManager()->hasFunctionSpaceOfType<TargetFunctionSpaceType3>(meshName))
      {
        // get target geometry field variable
        std::shared_ptr<TargetFieldVariableType3> geometryFieldTarget = std::make_shared<TargetFieldVariableType3>(
          this->context_.meshManager()->functionSpace<TargetFunctionSpaceType3>(meshName)->geometryField());

        LOG(DEBUG) << "transfer geometry field to linear mesh, " << geometryFieldSource->functionSpace()->meshName() << " -> "
          << geometryFieldTarget->functionSpace()->meshName();
        LOG(DEBUG) << StringUtility::demangle(typeid(SourceFunctionSpaceType).name()) << " -> " << StringUtility::demangle(typeid(TargetFunctionSpaceType3).name());

        // perform the mapping
        DihuContext::mappingBetweenMeshesManager()->template prepareMapping<SourceFieldVariableType,TargetFieldVariableType3>(geometryFieldSource, geometryFieldTarget, -1);

        // map the whole geometry field (all components), do not avoid copy
        DihuContext::mappingBetweenMeshesManager()->template map<SourceFieldVariableType,TargetFieldVariableType3>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
        DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<SourceFieldVariableType,TargetFieldVariableType3>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
      }

      // for second order composite meshes
      using TargetFunctionSpaceType4 = ::FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>;
      using TargetFieldVariableType4 = FieldVariable::FieldVariable<TargetFunctionSpaceType4,3>;

      // if the mesh name corresponds to a quadratic mesh
      if (this->context_.meshManager()->hasFunctionSpaceOfType<TargetFunctionSpaceType4>(meshName))
      {
        // get target geometry field variable
        std::shared_ptr<TargetFieldVariableType4> geometryFieldTarget = std::make_shared<TargetFieldVariableType4>(
          this->context_.meshManager()->functionSpace<TargetFunctionSpaceType4>(meshName)->geometryField());

        LOG(DEBUG) << "transfer geometry field to quadratic mesh, " << geometryFieldSource->functionSpace()->meshName() << " -> "
          << geometryFieldTarget->functionSpace()->meshName();
        LOG(DEBUG) << StringUtility::demangle(typeid(SourceFunctionSpaceType).name()) << " -> " << StringUtility::demangle(typeid(TargetFunctionSpaceType4).name());

        // perform the mapping
        DihuContext::mappingBetweenMeshesManager()->template prepareMapping<SourceFieldVariableType,TargetFieldVariableType4>(geometryFieldSource, geometryFieldTarget, -1);

        // map the whole geometry field (all components), do not avoid copy
        DihuContext::mappingBetweenMeshesManager()->template map<SourceFieldVariableType,TargetFieldVariableType4>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
        DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<SourceFieldVariableType,TargetFieldVariableType4>(geometryFieldSource, geometryFieldTarget, -1, -1, false);
      }
    }
  }

  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_+std::string("_map_geometry"));
}


//! get a reference to the DynamicHyperelasticitySolverType
template<typename MeshType,typename Term,bool withLargeOutputFiles>
std::shared_ptr<typename MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::DynamicHyperelasticitySolverType> MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
dynamicHyperelasticitySolver()
{
  return dynamicHyperelasticitySolver_;
}

template<typename MeshType,typename Term,bool withLargeOutputFiles>
typename MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::Data &MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
data()
{
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the slot_connector_data_transfer class
template<typename MeshType,typename Term,bool withLargeOutputFiles>
std::shared_ptr<typename MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::SlotConnectorDataType> MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
getSlotConnectorData()
{
  return data_.getSlotConnectorData();
}

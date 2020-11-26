#include "specialized_solver/solid_mechanics/hyperelasticity/00_initialize.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"
#include "partition/mesh_partition/01_mesh_partition_structured.h"
#include "specialized_solver/solid_mechanics/hyperelasticity/02_petsc_callbacks.h"

namespace SpatialDiscretization
{

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
HyperelasticityInitialize(DihuContext context, std::string settingsKey) :
  context_(context[settingsKey]), data_(context_), pressureDataCopy_(context_), initialized_(false),
  endTime_(0), lastNorm_(0), secondLastNorm_(0), currentLoadFactor_(1.0), lastSolveSucceeded_(true), nNonZerosJacobian_(0)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }

  // parse options concerning jacobian
  useAnalyticJacobian_  = this->specificSettings_.getOptionBool("useAnalyticJacobian", true);
  useNumericJacobian_   = this->specificSettings_.getOptionBool("useNumericJacobian", true);
  nNonlinearSolveCalls_ = this->specificSettings_.getOptionInt("nNonlinearSolveCalls", 1, PythonUtility::Positive);
  loadFactorGiveUpThreshold_ = this->specificSettings_.getOptionDouble("loadFactorGiveUpThreshold", 1e-5, PythonUtility::Positive);

  // parse constant body force, a value of "None" yields the default value, (0,0,0)
  constantBodyForce_ = this->specificSettings_.template getOptionArray<double,3>("constantBodyForce", Vec3{0.0,0.0,0.0});

  if (!useAnalyticJacobian_ && !useNumericJacobian_)
  {
    LOG(WARNING) << "Cannot set both \"useAnalyticJacobian\" and \"useNumericJacobian\" to False, now using numeric jacobian.";
    useNumericJacobian_ = true;
  }

  // parse material parameters
  specificSettings_.getOptionVector("materialParameters", materialParameters_);

  LOG(DEBUG) << "HyperelasticitySolver: parsed parameters " << materialParameters_;

  if (materialParameters_.size() != Term::nMaterialParameters)
  {
    LOG(FATAL) << "A wrong number of material parameters was specified. Specified parameters: " << materialParameters_.size()
      << " (" << materialParameters_ << "), needed parameters by " << StringUtility::demangle(typeid(Term).name()) << ": "
      << Term::nMaterialParameters;
  }

  // initialize material parameters
  // set all PARAM(i) values to the values given by materialParameters
  SEMT::set_parameters<Term::nMaterialParameters>::to(materialParameters_);

  displacementsScalingFactor_ = specificSettings_.getOptionDouble("displacementsScalingFactor", 1.0);
  dumpDenseMatlabVariables_ = specificSettings_.getOptionBool("dumpDenseMatlabVariables", false);

  // for the dynamic equation
  if (nDisplacementComponents == 6)
  {
    density_                 = specificSettings_.getOptionDouble("density", 1.0, PythonUtility::Positive);
    timeStepWidth_           = specificSettings_.getOptionDouble("timeStepWidth", 1.0, PythonUtility::Positive);
    extrapolateInitialGuess_ = specificSettings_.getOptionBool("extrapolateInitialGuess", true);
  }

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
  this->outputWriterManagerPressure_.initialize(this->context_["pressure"], this->context_["pressure"].getPythonConfig());
  this->outputWriterManagerLoadIncrements_.initialize(this->context_["LoadIncrements"], this->context_["LoadIncrements"].getPythonConfig());
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
initialize()
{
  if (this->initialized_)
    return;

  LOG(DEBUG) << "initialize HyperelasticitySolver";
  assert(this->specificSettings_.pyObject());

  // create function space / mesh, the geometry is from the settings
  displacementsFunctionSpace_ = context_.meshManager()->functionSpace<DisplacementsFunctionSpace>(specificSettings_);

  pressureFunctionSpace_ = PressureFunctionSpaceCreator<typename PressureFunctionSpace::Mesh>::createPressureFunctionSpace(context_.meshManager(), displacementsFunctionSpace_);

  // initialize the data object
  // store mesh in data
  LOG(DEBUG) << "call setDisplacementsFunctionSpace of data";
  data_.setDisplacementsFunctionSpace(displacementsFunctionSpace_);
  data_.setPressureFunctionSpace(pressureFunctionSpace_);

  data_.initialize();
  pressureDataCopy_.initialize(data_.pressure(), data_.displacementsLinearMesh(), data_.velocitiesLinearMesh());
  pressureDataCopy_.setFunctionSpace(pressureFunctionSpace_);

  // create nonlinear solver PETSc context (snes)
  nonlinearSolver_ = this->context_.solverManager()->template solver<Solver::Nonlinear>(
    this->specificSettings_, this->displacementsFunctionSpace_->meshPartition()->mpiCommunicator());

  // initialize Dirichlet boundary conditions
  if (dirichletBoundaryConditions_ == nullptr)
  {
    dirichletBoundaryConditions_ = std::make_shared<DirichletBoundaryConditions<DisplacementsFunctionSpace,nDisplacementComponents>>(this->context_);
    dirichletBoundaryConditions_->initialize(this->specificSettings_, this->data_.functionSpace(), "dirichletBoundaryConditions");
  }

  // initialize Neumann boundary conditions
  if (neumannBoundaryConditions_ == nullptr)
  {
    typedef Quadrature::Gauss<3> QuadratureType;
    neumannBoundaryConditions_ = std::make_shared<NeumannBoundaryConditions<DisplacementsFunctionSpace,QuadratureType,3>>(this->context_);
    neumannBoundaryConditions_->initialize(this->specificSettings_, this->data_.functionSpace(), "neumannBoundaryConditions");
    neumannBoundaryConditions_->setDeformationGradientField(this->data_.deformationGradient());
  }

  // initialize fiber direction field
  if (Term::usesFiberDirection)
  {
    initializeFiberDirections();
  }

  // setup Petsc variables
  LOG(DEBUG) << "initialize Petsc Variables";
  initializePetscVariables();

  // parse load factors
  this->specificSettings_.getOptionVector("loadFactors", loadFactors_);

  if (!loadFactors_.empty())
  {
    if (fabs(loadFactors_.back() - 1.0) > 1e-12)
    {
      LOG(WARNING) << this->specificSettings_ << "[\"loadFactors\"]: Last load factor " << loadFactors_.back() << " is not 1.0.";
    }
  }

  // prepare load factors
  if (loadFactors_.empty())
    loadFactors_.push_back(1.0);

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("HyperelasticitySolver");

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
initializeFiberDirections()
{
  std::vector<std::string> fiberMeshNames;
  this->specificSettings_.template getOptionVector<std::string>("fiberMeshNames", fiberMeshNames);

  LOG(INFO) << "initializing fiber directions...";

  LOG(DEBUG) << "config: " << this->specificSettings_;
  LOG(DEBUG) << fiberMeshNames.size() << " fiberMeshNames: " << fiberMeshNames;

  // if no fiber meshes were specified, use the settings fiberDirection
  if (fiberMeshNames.empty())
  {
    Vec3 fiberDirection = this->specificSettings_.template getOptionArray<double,3>("fiberDirection", Vec3{0,0,0});

    if (fiberDirection[0] == 0 && fiberDirection[1] == 0 && fiberDirection[2] == 0)
    {
      // fiberDirection was not specified, check localFiberDirection
      Vec3 fiberDirectionInElement = this->specificSettings_.template getOptionArray<double,3>("fiberDirectionInElement", Vec3{0,0,1});

      const int nNodes1D = ::FunctionSpace::FunctionSpaceBaseDim<1,BasisFunction::LagrangeOfOrder<2>>::nNodesPerElement();
      assert(nNodes1D == 3);
      const int nDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();

      // loop over elements
      for (element_no_t elementNoLocal = 0; elementNoLocal < this->displacementsFunctionSpace_->nElementsLocal(); elementNoLocal++)
      {
        // get all node positions of nodes of this element
        std::array<Vec3,nDofsPerElement> geometry;
        this->displacementsFunctionSpace_->getElementGeometry(elementNoLocal, geometry);

        // get indices of element-local dofs
        std::array<dof_no_t,nDofsPerElement> dofNosLocal = this->displacementsFunctionSpace_->getElementDofNosLocal(elementNoLocal);

        // loop over nodes in element
        for (int k = 0; k < nNodes1D; k++)
        {
          for (int j = 0; j < nNodes1D; j++)
          {
            for (int i = 0; i < nNodes1D; i++)
            {
              int elementalNodeNo = k*nNodes1D*nNodes1D + j*nNodes1D + i;
              dof_no_t dofNoLocal = dofNosLocal[elementalNodeNo];

              Vec3 elementalX, elementalY, elementalZ;
              getElementalBasis(i, j, k, geometry, elementalX, elementalY, elementalZ);

              Vec3 fiberDirection = fiberDirectionInElement[0] * elementalX + fiberDirectionInElement[1] * elementalY + fiberDirectionInElement[2] * elementalZ;
              this->data_.fiberDirection()->setValue(dofNoLocal, fiberDirection, INSERT_VALUES);
            }
          }
        }
      }
    }
    else
    {
      // fiberDirection was specified
      LOG(DEBUG) << "displacements field variable type data: " << StringUtility::demangle(typeid(typename Data::DisplacementsFieldVariableType).name());
      LOG(DEBUG) << "displacements function space type: " << StringUtility::demangle(typeid(DisplacementsFunctionSpace).name());
      LOG(DEBUG) << "displacements function space mesh partition: " << *this->displacementsFunctionSpace_->meshPartition();

      std::vector<Vec3> fiberDirections(this->displacementsFunctionSpace_->nDofsLocalWithGhosts(), fiberDirection);
      this->data_.fiberDirection()->setValues(this->displacementsFunctionSpace_->meshPartition()->dofNosLocal(), fiberDirections);
    }
  }
  else
  {
    // loop over fiber mesh names
    for (int fiberNo = 0; fiberNo < fiberMeshNames.size(); fiberNo++)
    {
      // get fiber function space
      std::string fiberMeshName = fiberMeshNames[fiberNo];
      LOG(DEBUG) << "fiber " << fiberNo << "/" << fiberMeshNames.size() << ", mesh \"" << fiberMeshName << "\".";

      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = context_.meshManager()->functionSpace<FiberFunctionSpace>(fiberMeshName);

      LOG(DEBUG) << "create mapping  1D fiber -> 3D: \"" << fiberFunctionSpace->meshName() << "\" -> \"" << this->displacementsFunctionSpace_->meshName() << "\".";

      // initialize mapping 1D fiber -> 3D
      DihuContext::mappingBetweenMeshesManager()->template mappingBetweenMeshes<FiberFunctionSpace,DisplacementsFunctionSpace>(fiberFunctionSpace, this->displacementsFunctionSpace_);
    }

    using SourceFunctionSpaceType = FieldVariable::FieldVariable<FiberFunctionSpace,3>;
    using TargetFunctionSpaceType = DisplacementsFieldVariableType;

    // prepare the target mesh for the mapping, set all factors to zero
    std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,3>> direction;
    DihuContext::mappingBetweenMeshesManager()->template prepareMapping<SourceFunctionSpaceType,TargetFunctionSpaceType>(direction, this->data_.fiberDirection(), -1);

    // loop over fiber mesh names
    for (int fiberNo = 0; fiberNo < fiberMeshNames.size(); fiberNo++)
    {
      std::string fiberMeshName = fiberMeshNames[fiberNo];
      LOG(DEBUG) << "mesh \"" << fiberMeshName << "\".";

      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = context_.meshManager()->functionSpace<FiberFunctionSpace>(fiberMeshName);

      // define direction field variable on the fiber that will store the direction of the fiber
      std::vector<std::string> components({"x","y","z"});
      direction = fiberFunctionSpace->createFieldVariable<3>("direction", components);

      // set values of direction field variable
      int nDofsLocalWithoutGhosts = fiberFunctionSpace->nDofsLocalWithoutGhosts();
      std::vector<Vec3> geometryFieldValues;
      fiberFunctionSpace->geometryField().getValuesWithoutGhosts(geometryFieldValues);

      std::vector<Vec3> directionValues(nDofsLocalWithoutGhosts);

      // loop over local nodes
      for (dof_no_t dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts; dofNoLocal++)
      {
        dof_no_t index0 = std::max((dof_no_t)(0), dofNoLocal-1);
        dof_no_t index1 = std::min((dof_no_t)(nDofsLocalWithoutGhosts-1), (dof_no_t)(dofNoLocal+1));

        // get direction of 1D fiber
        Vec3 fiberDirection = -geometryFieldValues[index0] + geometryFieldValues[index1];
        directionValues[dofNoLocal] = fiberDirection;
      }

      direction->setValuesWithoutGhosts(directionValues);


      // transfer direction values from 1D fibers to 3D field, -1 means all components
      DihuContext::mappingBetweenMeshesManager()->template map<SourceFunctionSpaceType,TargetFunctionSpaceType>(
        direction, this->data_.fiberDirection(), -1, -1, false);
    }

    // finalize the mapping to the target mesh, compute final values by dividing by the factors
    DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<SourceFunctionSpaceType,TargetFunctionSpaceType>(direction, this->data_.fiberDirection(), -1, -1, false);
  }

  LOG(DEBUG) << "normalize fiber direction";

  auto tStart = std::chrono::steady_clock::now();

  // normalize entries
  std::vector<Vec3> valuesLocalWithoutGhosts;
  this->data_.fiberDirection()->getValuesWithoutGhosts(valuesLocalWithoutGhosts);

  for (dof_no_t dofNoLocal = 0; dofNoLocal < this->data_.fiberDirection()->nDofsLocalWithoutGhosts(); dofNoLocal++)
  {
    if (MathUtility::normSquared<3>(valuesLocalWithoutGhosts[dofNoLocal]) < 1e-10)
    {
      valuesLocalWithoutGhosts[dofNoLocal] = Vec3({0, 0, 1});
      LOG(DEBUG) << "dof " << dofNoLocal << ", fiberDirection was (0,0,0), set to (0,0,1)";
    }
    else
    {
      MathUtility::normalize<3>(valuesLocalWithoutGhosts[dofNoLocal]);
#ifndef NDEBUG
      LOG(DEBUG) << "dof " << dofNoLocal << ", fiberDirection normalized: " << valuesLocalWithoutGhosts[dofNoLocal]
        << " (norm: " << MathUtility::norm<3>(valuesLocalWithoutGhosts[dofNoLocal]) << ")";
#endif
    }
  }

  this->data_.fiberDirection()->setValuesWithoutGhosts(valuesLocalWithoutGhosts);
  this->data_.fiberDirection()->zeroGhostBuffer();

  this->data_.fiberDirection()->finishGhostManipulation();
  this->data_.fiberDirection()->startGhostManipulation();
  //this->data_.fiberDirection()->zeroGhostBuffer();
  //this->data_.fiberDirection()->finishGhostManipulation();

  auto tEnd = std::chrono::steady_clock::now();
  LOG_N_TIMES(1,INFO) << "done (" << std::chrono::duration_cast<std::chrono::milliseconds>(tEnd-tStart).count() << " ms)";
}


template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
getElementalBasis(int i, int j, int k,
                  const std::array<Vec3,27> &geometry,
                  Vec3 &elementalX, Vec3 &elementalY, Vec3 &elementalZ)
{
  const int nNodes1D = ::FunctionSpace::FunctionSpaceBaseDim<1,BasisFunction::LagrangeOfOrder<2>>::nNodesPerElement(); // 3

  // get node position of node (i,j,k)
  int elementalNodeNo = k*nNodes1D*nNodes1D + j*nNodes1D + i;

  Vec3 nodePosition = geometry[elementalNodeNo];

  // get vector between neighbouring nodes in x direction
  int elementalNodeNoNextX;
  int sign;
  if (i == nNodes1D-1)
  {
    elementalNodeNoNextX = k*nNodes1D*nNodes1D + j*nNodes1D + i-1;
    sign = -1;
  }
  else
  {
    elementalNodeNoNextX = k*nNodes1D*nNodes1D + j*nNodes1D + i+1;
    sign = 1;
  }

  Vec3 nodePositionNextX = geometry[elementalNodeNoNextX];
  elementalX = (-nodePosition + nodePositionNextX) * sign;

  // get vector between neighbouring nodes in y direction
  int elementalNodeNoNextY;
  if (j == nNodes1D-1)
  {
    elementalNodeNoNextY = k*nNodes1D*nNodes1D + (j-1)*nNodes1D + i;
    sign = -1;
  }
  else
  {
    elementalNodeNoNextY = k*nNodes1D*nNodes1D + (j+1)*nNodes1D + i;
    sign = 1;
  }

  Vec3 nodePositionNextY = geometry[elementalNodeNoNextY];
  elementalY = (-nodePosition + nodePositionNextY) * sign;

  // get vector between neighbouring nodes in z direction
  int elementalNodeNoNextZ;
  if (k == nNodes1D-1)
  {
    elementalNodeNoNextZ = (k-1)*nNodes1D*nNodes1D + j*nNodes1D + i;
    sign = -1;
  }
  else
  {
    elementalNodeNoNextZ = (k+1)*nNodes1D*nNodes1D + j*nNodes1D + i;
    sign = 1;
  }

  Vec3 nodePositionNextZ = geometry[elementalNodeNoNextZ];
  elementalZ = (-nodePosition + nodePositionNextZ) * sign;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::VecHyperelasticity>
HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
createPartitionedPetscVec(std::string name)
{
  LOG(DEBUG) << "createPartitionedPetscVec(" << name << ")";
  return std::make_shared<VecHyperelasticity>(
    displacementsFunctionSpace_->meshPartition(), pressureFunctionSpace_->meshPartition(), dirichletBoundaryConditions_, name);
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::MatHyperelasticity>
HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
createPartitionedPetscMat(std::string name)
{
  // determine number of non zero entries in matrix
  if (nNonZerosJacobian_ == 0)
    nNonZerosJacobian_ = materialDetermineNumberNonzerosInJacobian();

  int nNonZerosOffdiagonal = (int)nNonZerosJacobian_;
  int nNonZerosDiagonal = (int) nNonZerosJacobian_;

  //::Data::FiniteElementsBase<DisplacementsFunctionSpace,1>::getPetscMemoryParameters(nNonZerosDiagonal, nNonZerosOffdiagonal);

  //nNonZerosDiagonal = 100*MathUtility::sqr(nNonZerosDiagonal);
  //nNonZerosOffdiagonal = 5*MathUtility::sqr(nNonZerosOffdiagonal);

  LOG(INFO) << "Preallocation for matrix \"" << name << "\": diagonal nz: " << nNonZerosDiagonal << ", offdiagonal nz: " << nNonZerosOffdiagonal;

  return std::make_shared<MatHyperelasticity>(
    combinedVecSolution_, nNonZerosDiagonal, nNonZerosOffdiagonal, name);
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
unsigned int HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
materialDetermineNumberNonzerosInJacobian()
{
  unsigned int nNonZeros = 0;

  // get pointer to function space
  std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace = this->data_.displacementsFunctionSpace();
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace = this->data_.pressureFunctionSpace();

  const int D = 3;  // dimension
  const int nDisplacementsDofsPerElement = DisplacementsFunctionSpace::nDofsPerElement();
  const int nPressureDofsPerElement = PressureFunctionSpace::nDofsPerElement();
  const int nElementsLocal = displacementsFunctionSpace->nElementsLocal();

  nNonZeros = nElementsLocal * MathUtility::sqr(nDisplacementsDofsPerElement * D);

  if (nDisplacementComponents == 6)
  {
    nNonZeros += nElementsLocal * MathUtility::sqr(nDisplacementsDofsPerElement * D) * 3;
  }

  if (Term::isIncompressible)
  {
    nNonZeros += nElementsLocal * nPressureDofsPerElement * nDisplacementsDofsPerElement * D * 2;
    nNonZeros += nElementsLocal * nPressureDofsPerElement;
  }

  return nNonZeros;
}

//! get the precomputed external virtual work
template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
Vec HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
externalVirtualWork()
{
  return externalVirtualWorkDead_;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
initializePetscVariables()
{
  /*
   * jacobian matrix layout (for one process and nDisplacementComponents==3):
   *  (U U U P)
   *  (U U U P)
   *  (U U U P)
   *  (P P P 0)
   */
  // prepare for data structures without ghost dofs, these are normal Petsc Mat's and Vec's for which all solvers are available

  // solution vector
  combinedVecSolution_ = createPartitionedPetscVec("combinedSolution");

  // residual vector
  combinedVecResidual_ = createPartitionedPetscVec("combinedResidual");

  // vector for the external virtual work contribution that does not depend on u, δW_ext,dead (this is the same as δW_ext for static case)
  combinedVecExternalVirtualWorkDead_ = createPartitionedPetscVec("combinedVecExternalVirtualWorkDead");

  // get number of local entries of the vectors, this will be the number of rows and columns of the matrix
  int nMatrixRowsLocal = combinedVecSolution_->nEntriesLocal();

  // output
  LOG(DEBUG) << "n dofs displacements: " << nDisplacementComponents << "*" << displacementsFunctionSpace_->nDofsLocalWithoutGhosts()
    << " = " << nDisplacementComponents*displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
  LOG(DEBUG) << "n dofs pressure: " << pressureFunctionSpace_->nDofsLocalWithoutGhosts() << ", total: "
    << displacementsFunctionSpace_->nDofsLocalWithoutGhosts() * nDisplacementComponents + pressureFunctionSpace_->nDofsLocalWithoutGhosts();
  LOG(DEBUG) << "n BC dofs: ";

  LOG(DEBUG) << "  ux: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[0].dofNosLocal.size();
  LOG(DEBUG) << "  uy: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[1].dofNosLocal.size();
  LOG(DEBUG) << "  uz: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[2].dofNosLocal.size();

  if (nDisplacementComponents == 6)
  {
    LOG(DEBUG) << "  vx: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[3].dofNosLocal.size();
    LOG(DEBUG) << "  vy: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[4].dofNosLocal.size();
    LOG(DEBUG) << "  vz: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[5].dofNosLocal.size();
  }
  LOG(DEBUG) << "number of non-BC dofs total: " << nMatrixRowsLocal;

  // create matrix with same dof mapping as vectors
  //std::shared_ptr<::FunctionSpace::Generic> genericFunctionSpace = context_.meshManager()->createGenericFunctionSpace(nMatrixRowsLocal, displacementsFunctionSpace_->meshPartition(), "genericMesh");

  combinedMatrixJacobian_ = createPartitionedPetscMat("combinedJacobian");

  solverMatrixAdditionalNumericJacobian_ = PETSC_NULL;

  // if both numeric and analytic jacobian are used, create additional matrix that will hold the numeric jacobian
  if (useNumericJacobian_ && useAnalyticJacobian_)
  {
    combinedMatrixAdditionalNumericJacobian_ = createPartitionedPetscMat("combinedJacobianNumeric");

    solverMatrixAdditionalNumericJacobian_ = combinedMatrixAdditionalNumericJacobian_->valuesGlobal();
  }

  // extract the Petsc Vec's of the PartitionedPetscVecForHyperelasticity objects
  LOG(DEBUG) << "get the internal vectors";
  solverMatrixJacobian_ = combinedMatrixJacobian_->valuesGlobal();
  solverVariableSolution_ = combinedVecSolution_->valuesGlobal();
  solverVariableResidual_ = combinedVecResidual_->valuesGlobal();
  externalVirtualWorkDead_ = combinedVecExternalVirtualWorkDead_->valuesGlobal();

  // create vector with all zeros in it, this is needed for zeroing the diagonal of the stiffness matrix for initialization in evaluateAnalyticJacobian
  PetscErrorCode ierr;
  ierr = VecDuplicate(solverVariableResidual_, &zeros_); CHKERRV(ierr);
  ierr = VecZeroEntries(zeros_); CHKERRV(ierr);
  ierr = VecDuplicate(solverVariableResidual_, &lastSolution_); CHKERRV(ierr);
  ierr = VecDuplicate(solverVariableResidual_, &bestSolution_); CHKERRV(ierr);

  LOG(DEBUG) << "for debugging: " << combinedVecSolution_->getString();

  // compute the external virtual work, because it is constant throughout the solution process
  materialComputeExternalVirtualWorkDead();

  // output pointer values for debugging
  LOG(DEBUG) << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.";
  LOG(DEBUG) << "pointer value solverVariableResidual_: " << solverVariableResidual_;
  LOG(DEBUG) << "pointer value solverVariableSolution_: " << solverVariableSolution_;
  LOG(DEBUG) << "pointer value solverMatrixJacobian_:   " << solverMatrixJacobian_;
  LOG(DEBUG) << "pointer value solverMatrixAdditionalNumericJacobian_:  " << solverMatrixAdditionalNumericJacobian_;
  LOG(DEBUG) << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.";

  // if a numeric jacobian is in use, assemble the matrix
  // depending on if also the analytic jacobian will be computed, the storage for the numeric jacobian is either combinedMatrixJacobian_ or combinedMatrixAdditionalNumericJacobian_
  if (useNumericJacobian_)
  {
    if (useAnalyticJacobian_)
    {
      // assemble matrix, new nonzeros are allocated later
      combinedMatrixAdditionalNumericJacobian_->assembly(MAT_FINAL_ASSEMBLY);
    }
    else
    {
      // assemble matrix, new nonzeros are allocated later
      combinedMatrixJacobian_->assembly(MAT_FINAL_ASSEMBLY);
    }
  }

  if (useAnalyticJacobian_)
  {
    // jacobian matrix is already preallocated, but there might be even more entries required
    ierr = MatSetOption(solverMatrixJacobian_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRV(ierr);

    // assemble matrix and nonzeros structure
    evaluateAnalyticJacobian(solverVariableSolution_, solverMatrixJacobian_);

    // print info about preallocation, this can also be done by -mat_view ::ascii_info_detail
    //ierr = MatView(solverMatrixJacobian_, PETSC_VIEWER_ASCII_INFO_DETAIL); CHKERRV(ierr);

    // output the jacobian matrix for debugging
    LOG(DEBUG) << "initial analytic jacobian matrix: ";
    dumpJacobianMatrix(solverMatrixJacobian_);
  }

  // assign all callback functions
  this->initializePetscCallbackFunctions();

  // set solution vector to zero or initial value
  this->initializeSolutionVariable();
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
initializeSolutionVariable()
{
  // set variable to all zero and dirichlet boundary condition value
  LOG(DEBUG) << "zeroEntries, representation: " << combinedVecSolution_->currentRepresentation();
  combinedVecSolution_->zeroEntries();


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
  }
  else
  {
    // input is already only the local dofs, use all
    const int nDofsLocal = displacementsFunctionSpace->nDofsLocalWithoutGhosts();
    this->specificSettings_.template getOptionVector<Vec3>("initialValuesDisplacements", nDofsLocal, localValuesDisplacements);
  }
  VLOG(1) << "set initial values for displacements to " << localValuesDisplacements;
  VLOG(1) << "set initial values for velocities to " << localValuesVelocities;

  // set the first component of the solution variable by the given values
  this->data_.displacements()->setValuesWithoutGhosts(localValuesDisplacements);

  // set displacement entries in combinedVecSolution_
  int nDofsLocalWithoutGhosts = displacementsFunctionSpace->nDofsLocalWithoutGhosts();
  std::vector<double> localValues(nDofsLocalWithoutGhosts);

  combinedVecSolution_->startGhostManipulation();

  // set displacement entries in combinedVecSolution_
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    for (int entryNo = 0; entryNo < nDofsLocalWithoutGhosts; entryNo++)
    {
      localValues[entryNo] = localValuesDisplacements[entryNo][componentNo];
    }

    combinedVecSolution_->setValues(componentNo, nDofsLocalWithoutGhosts, displacementsFunctionSpace->meshPartition()->dofNosLocal().data(), localValues.data());
  }

  // assemble vector
  combinedVecSolution_->zeroGhostBuffer();
  combinedVecSolution_->finishGhostManipulation();

  LOG(DEBUG) << "values: " << PetscUtility::getStringVector(solverVariableSolution_);
  LOG(DEBUG) << "after initialization: " << combinedVecSolution_->getString();
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
setDisplacementsAndPressureFromCombinedVec(Vec x, std::shared_ptr<DisplacementsFieldVariableType> u,
                                           std::shared_ptr<PressureFieldVariableType> p)
{
  assert(nDisplacementComponents == 3);

  // copy entries of combined vector x to this->data_.displacements() and this->data_.pressure()
  std::vector<double> values;

  if (VLOG_IS_ON(1))
  {
    PetscUtility::getVectorEntries(x, values);
    VLOG(1) << "setDisplacementsAndPressureFromCombinedVec, x=" << PetscUtility::getStringVector(x);
  }

  bool backupVecs = false;
  if (x != solverVariableSolution_)
  {
    backupVecs = true;
  }

  // move the values of x to variable solverVariableSolution_
  if (backupVecs)
  {
    // this happens in computation of the numeric jacobian
    VecSwap(x, solverVariableSolution_);
  }

  if (!u && !p)
  {
    u = this->data_.displacements();
    if (Term::isIncompressible)
      p = this->data_.pressure();   // p is only needed for incompressible formulation
  }

  // set displacement entries
  u->zeroGhostBuffer();
  if (p)
  {
    p->zeroGhostBuffer();
  }
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    int nEntries = displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
    values.resize(nEntries);
    combinedVecSolution_->getValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());

    u->setValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
  }
  u->finishGhostManipulation();
  u->startGhostManipulation();

  // set pressure entries
  if (p)
  {
    int nEntries = pressureFunctionSpace_->nDofsLocalWithoutGhosts();
    values.resize(nEntries);

    const int pressureComponent = nDisplacementComponents;
    combinedVecSolution_->getValues(pressureComponent, nEntries, pressureFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    p->setValues(0, nEntries, pressureFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    p->finishGhostManipulation();
    p->startGhostManipulation();
  }

  // undo the backup operation
  if (backupVecs)
  {
    VecSwap(x, solverVariableSolution_);
  }

  //VLOG(1) << *u;
  //VLOG(1) << *p;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
setDisplacementsVelocitiesAndPressureFromCombinedVec(Vec x,
                                                     std::shared_ptr<DisplacementsFieldVariableType> u,
                                                     std::shared_ptr<DisplacementsFieldVariableType> v,
                                                     std::shared_ptr<PressureFieldVariableType> p)
{
  assert(nDisplacementComponents == 6);

  // copy entries of combined vector x to this->data_.displacements() and this->data_.pressure()
  std::vector<double> values;

  if (VLOG_IS_ON(1))
  {
    PetscUtility::getVectorEntries(x, values);
    VLOG(1) << "setDisplacementsVelocitiesAndPressureFromCombinedVec, x=" << PetscUtility::getStringVector(x);
  }

  bool backupVecs = false;
  if (x != solverVariableSolution_)
  {
    backupVecs = true;
  }

  // move the values of x to variable solverVariableSolution_
  if (backupVecs)
  {
    // this happens in computation of the numeric jacobian
    VecSwap(x, solverVariableSolution_);
  }

  if (!u && !v && !p)
  {
    u = this->data_.displacements();
    v = this->data_.velocities();

    if (Term::isIncompressible)
      p = this->data_.pressure();   // p is only needed for incompressible formulation
  }

  // set displacement entries
  u->zeroGhostBuffer();
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    int nEntries = displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
    values.resize(nEntries);
    combinedVecSolution_->getValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());

    if (VLOG_IS_ON(1))
      LOG(DEBUG) << "setDisplacementsVelocitiesAndPressureFromCombinedVec, " << nEntries << " u values: " << values;

    u->setValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
  }
  u->finishGhostManipulation();
  u->startGhostManipulation();

  // set velocity entries
  if (v)
  {
    v->zeroGhostBuffer();
    for (int componentNo = 0; componentNo < 3; componentNo++)
    {
      int nEntries = displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
      values.resize(nEntries);
      combinedVecSolution_->getValues(3 + componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());

    if (VLOG_IS_ON(1))
      VLOG(1) << "setDisplacementsVelocitiesAndPressureFromCombinedVec, " << nEntries << " v values: " << values;

      v->setValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    }
    v->finishGhostManipulation();
    v->startGhostManipulation();
  }
  // set pressure entries
  if (p)
  {
    p->zeroGhostBuffer();
    int nEntries = pressureFunctionSpace_->nDofsLocalWithoutGhosts();
    values.resize(nEntries);

    const int pressureComponent = nDisplacementComponents;
    combinedVecSolution_->getValues(pressureComponent, nEntries, pressureFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());

    if (VLOG_IS_ON(1))
      VLOG(1) << "setDisplacementsVelocitiesAndPressureFromCombinedVec, " << nEntries << " p values: " << values;

    p->setValues(0, nEntries, pressureFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    p->finishGhostManipulation();
    p->startGhostManipulation();
  }

  // undo the backup operation
  if (backupVecs)
  {
    VecSwap(x, solverVariableSolution_);
  }

  //VLOG(1) << *u;
  //VLOG(1) << *p;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
setUVP(Vec x)
{
  if (nDisplacementComponents == 3)
  {
    setDisplacementsAndPressureFromCombinedVec(x);
  }
  else if (nDisplacementComponents == 6)
  {
    setDisplacementsVelocitiesAndPressureFromCombinedVec(x);
  }
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
dumpJacobianMatrix(Mat jac)
{
  if (!dumpDenseMatlabVariables_)
    return;

  static int evaluationNo = 0;  // counter how often this function was called

  std::stringstream filename;
  filename << "out/jac_" << std::setw(3) << std::setfill('0') << evaluationNo << "r" <<  DihuContext::nRanksCommWorld();

  evaluationNo++;

  // if there are both numeric and analytic jacobian in use, this is the numeric jacobian
  if (jac == solverMatrixAdditionalNumericJacobian_)
  {
    filename << "_numeric";
    combinedMatrixAdditionalNumericJacobian_->dumpMatrixGlobalNatural(filename.str());
  }
  else if (jac == solverMatrixJacobian_)
  {
    // this is the normal jacobian, either numeric or analytic, if only one of both is in use

    // if both matrices are used
    if (solverMatrixAdditionalNumericJacobian_ != PETSC_NULL)
    {
      filename << "_analytic";
      combinedMatrixJacobian_->dumpMatrixGlobalNatural(filename.str());

      evaluationNo--;

      // if matrix is not all zeros
      double numericJacobianNorm = 0;
      MatNorm(solverMatrixAdditionalNumericJacobian_, NORM_1, &numericJacobianNorm);

      LOG(DEBUG) << "numericJacobianNorm: " << numericJacobianNorm;

      if (numericJacobianNorm > 1e-5)
      {
        // compute difference between analytic and numeric jacobian
        Mat difference;
        PetscErrorCode ierr;
        ierr = MatDuplicate(solverMatrixAdditionalNumericJacobian_, MAT_COPY_VALUES, &difference); CHKERRV(ierr);
        ierr = MatCopy(solverMatrixAdditionalNumericJacobian_, difference, SAME_NONZERO_PATTERN); CHKERRV(ierr);
        MatAXPY(difference, -1, solverMatrixJacobian_, DIFFERENT_NONZERO_PATTERN);

        double norm1 = 0;
        double normF = 0;
        double normInf = 0;
        MatNorm(difference, NORM_1, &norm1);
        MatNorm(difference, NORM_FROBENIUS, &normF);
        MatNorm(difference, NORM_INFINITY, &normInf);
        LOG(INFO) << "difference between analytic and numeric jacobian matrices: "
          << "1-norm: " << norm1 << ", frobenius norm: " << normF << ", infinity norm: " << normInf;

        // compute differences for submatrices
        int nRows = 2;
        if (nDisplacementComponents == 6)
          nRows = 3;
        if (!Term::isIncompressible)    // compressible formulation does not have pressure component
          nRows--;

        for (int i = 0; i < nRows; i++)
        {
          for (int j = 0; j < nRows; j++)
          {
            Mat analyticJacobianSubmatrix = combinedMatrixJacobian_->getSubmatrix(i,j);
            Mat numericJacobianSubmatrix = combinedMatrixAdditionalNumericJacobian_->getSubmatrix(i,j);

            Mat differenceSubmatrix;
            ierr = MatDuplicate(analyticJacobianSubmatrix, MAT_COPY_VALUES, &differenceSubmatrix); CHKERRV(ierr);
            ierr = MatCopy(analyticJacobianSubmatrix, differenceSubmatrix, SAME_NONZERO_PATTERN); CHKERRV(ierr);
            MatAXPY(differenceSubmatrix, -1, numericJacobianSubmatrix, DIFFERENT_NONZERO_PATTERN);

            MatNorm(differenceSubmatrix, NORM_1, &norm1);
            MatNorm(differenceSubmatrix, NORM_FROBENIUS, &normF);
            MatNorm(differenceSubmatrix, NORM_INFINITY, &normInf);
            LOG(INFO) << " submatrix (" << i << "," << j << "): 1-norm: " << norm1 << ", frobenius norm: " << normF << ", infinity norm: " << normInf;
          }
        }

        if (norm1 > 1)
          LOG(ERROR) << "norm mismatch";
      }
    }
    else
    {
      combinedMatrixJacobian_->dumpMatrixGlobalNatural(filename.str());
    }
  }
  else
  {
    LOG(ERROR) << "Could not output jacobian matrix " << jac << " (solverMatrixJacobian_: " << solverMatrixJacobian_ << ", solverMatrixAdditionalNumericJacobian_: " << solverMatrixAdditionalNumericJacobian_ << ")";
  }
}

//! get the PartitionedPetsVec for the residual and result of the nonlinear function
template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::VecHyperelasticity> HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
combinedVecResidual()
{
  return this->combinedVecResidual_;
}

//! get the PartitionedPetsVec for the solution
template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::VecHyperelasticity> HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
combinedVecSolution()
{
  return this->combinedVecSolution_;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
reset()
{
  this->initialized_ = false;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
typename HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::Data &HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
data()
{
  return data_;
}

//! get a pointer to the dirichlet boundary conditions object
template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
std::shared_ptr<DirichletBoundaryConditions<typename HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::DisplacementsFunctionSpace,nDisplacementComponents>>
HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
dirichletBoundaryConditions()
{
  return dirichletBoundaryConditions_;
}

//! get a pointer to the neumann boundary conditions object
template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
std::shared_ptr<NeumannBoundaryConditions<typename HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::DisplacementsFunctionSpace,Quadrature::Gauss<3>,3>>
HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
neumannBoundaryConditions()
{
  return neumannBoundaryConditions_;
}

//! set new neumann bc's = traction for the next solve
template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
updateNeumannBoundaryConditions(std::shared_ptr<NeumannBoundaryConditions<typename HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::DisplacementsFunctionSpace,Quadrature::Gauss<3>,3>> newNeumannBoundaryConditions)
{
  neumannBoundaryConditions_ = newNeumannBoundaryConditions;

  // compute new value for the rhs, δW_ext,dead = int_Ω B^L * phi^L * phi^M * δu^M dx + int_∂Ω T^L * phi^L * phi^M * δu^M dS
  materialComputeExternalVirtualWorkDead();
}

//! set new dirichlet boundary condition values for existing dofs
template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
updateDirichletBoundaryConditions(std::vector<std::pair<global_no_t,std::array<double,3>>> newDirichletBCValues)
{
  bool inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
  combinedVecSolution_->updateDirichletBoundaryConditions(newDirichletBCValues, inputMeshIsGlobal);
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
addDirichletBoundaryConditions(std::vector<typename DirichletBoundaryConditions<DisplacementsFunctionSpace,nDisplacementComponents>::ElementWithNodes> &boundaryConditionElements, bool overwriteBcOnSameDof)
{
  LOG(DEBUG) << "addDirichletBoundaryConditions, Term: " << StringUtility::demangle(typeid(Term).name());
  if (!Term::isIncompressible)
  {
    LOG(DEBUG) << "addDirichletBoundaryConditions on compressible material";
  }

  dirichletBoundaryConditions_->addBoundaryConditions(boundaryConditionElements, overwriteBcOnSameDof);

  // an incompressible material has 3+1 (static problem) or 6+1 (dynamic problem) components (displacements+velocities+pressure)
  // a compressible material has only 3 (static problem) or 6 (dynamic problem) components (no pressure)

  const int nComponents = nDisplacementComponents + (Term::isIncompressible? 1 : 0);

  // save previous values
  std::array<std::vector<double>, nComponents> combinedSolutionValues;
  std::array<std::vector<double>, nComponents> combinedResidualValues;
  std::array<std::vector<double>, nComponents> combinedVecExternalVirtualWorkDeadValues;

  int nDisplacementsVelocityValues = this->displacementsFunctionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
  int nPressureValues = this->pressureFunctionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
  std::vector<int> indices(nDisplacementsVelocityValues);
  std::iota(indices.begin(), indices.end(), 0);

  // loop over displacement components (and velocity components, if dynamic) and one pressure component
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    int nValues = nDisplacementsVelocityValues;
    if (componentNo == nDisplacementComponents)
      nValues = nPressureValues;

    // solution vector
    combinedSolutionValues[componentNo].resize(nValues);
    combinedVecSolution_->getValues(componentNo, nValues, indices.data(), combinedSolutionValues[componentNo].data());

    // residual vector
    //combinedResidualValues[componentNo].resize(nValues);
    //combinedVecResidual_->getValues(componentNo, nValues, indices.data(), combinedResidualValues[componentNo].data());

    // vector for the external virtual work contribution that does not depend on u, δW_ext,dead (this is the same as δW_ext for static case)
    combinedVecExternalVirtualWorkDeadValues[componentNo].resize(nValues);
    combinedVecExternalVirtualWorkDead_->getValues(componentNo, nValues, indices.data(), combinedVecExternalVirtualWorkDeadValues[componentNo].data());
  }

  // remove generic meshes "genericMesh" and "genericMeshForMatrixcombinedJacobian"
  DihuContext::meshManager()->deleteFunctionSpace("genericMesh");
  DihuContext::meshManager()->deleteFunctionSpace("genericMeshForMatrixcombinedJacobian");

  // create new vectors and matrices with the updated boundary conditions
  initializePetscVariables();

  // restore previous values, except for the new Dirichlet values
  // loop over displacement components and one pressure component
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    int nValues = nDisplacementsVelocityValues;
    if (componentNo == nDisplacementComponents)
      nValues = nPressureValues;

    // solution vector
    combinedVecSolution_->setValues(componentNo, nValues, indices.data(), combinedSolutionValues[componentNo].data());

    // residual vector
    //combinedVecResidual_->setValues(componentNo, nValues, indices.data(), combinedResidualValues[componentNo].data());

    // vector for the external virtual work contribution that does not depend on u, δW_ext,dead (this is the same as δW_ext for static case)
    combinedVecExternalVirtualWorkDead_->setValues(componentNo, nValues, indices.data(), combinedVecExternalVirtualWorkDeadValues[componentNo].data());
  }

  combinedVecSolution_->startGhostManipulation();
  combinedVecSolution_->finishGhostManipulation();

  //combinedVecResidual_->startGhostManipulation();
  //combinedVecResidual_->finishGhostManipulation();

  //combinedVecExternalVirtualWorkDead_->startGhostManipulation();
  //combinedVecExternalVirtualWorkDead_->finishGhostManipulation();

  PetscErrorCode ierr;
  ierr = VecAssemblyBegin(solverVariableSolution_); CHKERRV(ierr);
  ierr = VecAssemblyEnd(solverVariableSolution_); CHKERRV(ierr);
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
Vec HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
currentState()
{
  return solverVariableSolution_;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
std::string HyperelasticityInitialize<Term,withLargeOutput,MeshType,nDisplacementComponents>::
getString(Vec x)
{
  if (x == solverVariableSolution_)
  {
    return combinedVecSolution_->getString();
  }
  else if (x == solverVariableResidual_)
  {
    return combinedVecResidual_->getString();
  }
  else if (x == externalVirtualWorkDead_)
  {
    return combinedVecExternalVirtualWorkDead_->getString();
  }
  else
  {
    LOG(FATAL) << "this should not be called";
  }
  return std::string("no getString representation");
}


} // namespace SpatialDiscretization

#include "spatial_discretization/finite_element_method/00_base.h"

#include <Python.h>  // this has to be the first included header
#include <iostream>
#include <petscksp.h>
#include <memory>
#include <cassert>

#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"

#include "mesh/structured_regular_fixed.h"
#include "basis_function/lagrange.h"
#include "mesh/mesh_manager/mesh_manager.h"
#include "solver/solver_manager.h"
#include "solver/linear.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

namespace SpatialDiscretization
{

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
FiniteElementMethodBase(DihuContext context, std::shared_ptr<FunctionSpaceType> functionSpace) :
  context_(context["FiniteElementMethod"]), data_(context_), specificSettings_(context_.getPythonConfig()), initialized_(false)
{
  LOG(DEBUG) << "FiniteElementMethodBase constructor, context: " << this->context_.getPythonConfig();
  outputWriterManager_.initialize(context_, specificSettings_);

  // Create mesh or retrieve existing mesh from meshManager. This already creates meshPartition in functionSpace.initialize(), see function_space/03_function_space_partition_structured.tpp
  if (!functionSpace)
  {
    LOG(DEBUG) << "FiniteElementMethodBase constructor, create new function space from settings";
    functionSpace = context_.meshManager()->functionSpace<FunctionSpaceType>(specificSettings_);
  }
  else
  {
    LOG(DEBUG) << "FiniteElementMethodBase constructor, use given functionSpace \"" << functionSpace->meshName() << "\"";
    if (VLOG_IS_ON(1))
    {
      VLOG(1) << "geometry field: " << functionSpace->geometryField();
    }
  }

  // store mesh in data
  data_.setFunctionSpace(functionSpace);
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
std::shared_ptr<FunctionSpaceType> FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
functionSpace()
{
  return data_.functionSpace();
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
Data::FiniteElements<FunctionSpaceType,nComponents,Term> &FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
data()
{
  return data_;
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
setRankSubset(Partition::RankSubset rankSubset)
{
  data_.setRankSubset(rankSubset);
}
 
template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
initialize()
{
  // do not initialize if it was called already
  if (initialized_)
    return;

  data_.initialize();

  if (specificSettings_.hasKey("updatePrescribedValuesFromSolution"))
  {
    updatePrescribedValuesFromSolution_ = specificSettings_.getOptionBool("updatePrescribedValuesFromSolution", false);
    LOG(DEBUG) << "set updatePrescribedValuesFromSolution = " << updatePrescribedValuesFromSolution_;
  }

  // assemble stiffness matrix
  Control::PerformanceMeasurement::start("durationSetStiffnessMatrix");

  // initialize spatial parameter prefactor
  prefactor_.initialize(specificSettings_, "prefactor", 1.0, this->data_.functionSpace());

  // compute the stiffness matrix
  setStiffnessMatrix();

  // save the stiffness matrix also in the other slot, that will not be overwritten by applyBoundaryConditions
  PetscErrorCode ierr = MatDuplicate(this->data_.stiffnessMatrix()->valuesGlobal(), MAT_COPY_VALUES,
                                     &this->data_.stiffnessMatrixWithoutBc()->valuesGlobal()); CHKERRV(ierr);
  this->data_.stiffnessMatrixWithoutBc()->assembly(MAT_FINAL_ASSEMBLY);

  Control::PerformanceMeasurement::stop("durationSetStiffnessMatrix");

  if (updatePrescribedValuesFromSolution_)
  {
    PetscUtility::dumpMatrix("stiffnessmatrix_w", "matlab", this->data_.stiffnessMatrixWithoutBc()->valuesGlobal(), MPI_COMM_WORLD);
    PetscUtility::dumpMatrix("stiffnessmatrix", "matlab", this->data_.stiffnessMatrix()->valuesGlobal(), MPI_COMM_WORLD);
  }

  // set the rhs
  Control::PerformanceMeasurement::start("durationSetRightHandSide");
  setRightHandSide();
  Control::PerformanceMeasurement::stop("durationSetRightHandSide");

  // apply boundary conditions
  Control::PerformanceMeasurement::start("durationAssembleBoundaryConditions");
  this->applyBoundaryConditions();
  Control::PerformanceMeasurement::stop("durationAssembleBoundaryConditions");

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("FiniteElementMethod");
  DihuContext::solverStructureVisualizer()->setSlotConnectorData(getSlotConnectorData());

  initialized_ = true;
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
reset()
{
  data_.reset();
  initialized_ = false;
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
run()
{
  initialize();
  applyBoundaryConditions();
  solve();
  data_.print();

  callOutputWriter(-1, 0.0);
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  outputWriterManager_.writeOutput(data_, timeStepNo, currentTime, callCountIncrement);
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
solve()
{
  // solve linear system k*d=f for d
  LOG(TRACE) << "FiniteElementMethod::solve";

  // if equation was set to none, do not solve the problem (this is for unit tests that don't test for solution)
  if (std::is_same<Term,Equation::None>::value)
  {
    // set solution to zero
    data_.solution()->zeroEntries();
    return;
  }

  // get stiffness matrix
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix = data_.stiffnessMatrix();

  // assemble matrix such that all entries are at their place
  stiffnessMatrix->assembly(MAT_FINAL_ASSEMBLY);
  
  // get linear solver context from solver manager
  std::shared_ptr<Solver::Linear> linearSolver = this->context_.solverManager()->template solver<Solver::Linear>(
    this->specificSettings_, this->data_.functionSpace()->meshPartition()->mpiCommunicator());
  std::shared_ptr<KSP> ksp = linearSolver->ksp();
  assert(ksp != nullptr);

  // set matrix used for linear system and preconditioner to ksp context
  PetscErrorCode ierr;
  ierr = KSPSetOperators(*ksp, stiffnessMatrix->valuesGlobal(), stiffnessMatrix->valuesGlobal()); CHKERRV(ierr);

  VLOG(1) << "rhs: " << *data_.rightHandSide();
  VLOG(1) << "stiffnessMatrix: " << *stiffnessMatrix;

  // initialize coordinates for PETSc geometric multi-grid solvers
  setInformationToPreconditioner();

  // non-zero initial values
#if 0  
  PetscScalar scalar = 0.5;
  ierr = VecSet(data_.solution()->values(), scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRV(ierr);
#endif

  LOG(DEBUG) << "solve...";

  // solve the system
  linearSolver->solve(data_.rightHandSide()->valuesGlobal(), data_.solution()->valuesGlobal(), "Solution obtained");

  data_.solution()->setRepresentationGlobal();
  data_.solution()->startGhostManipulation();
  data_.solution()->zeroGhostBuffer();
  data_.solution()->finishGhostManipulation();
  
  
  VLOG(1) << "solution: " << *data_.solution();
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
setInformationToPreconditioner()
{
  // get linear solver context from solver manager
  std::shared_ptr<Solver::Linear> linearSolver = this->context_.solverManager()->template solver<Solver::Linear>(
    this->specificSettings_, this->data_.functionSpace()->meshPartition()->mpiCommunicator());
  std::shared_ptr<KSP> ksp = linearSolver->ksp();
  assert(ksp != nullptr);
  PetscErrorCode ierr;

  // get preconditioner
  PC pc;
  ierr = KSPGetPC(*ksp, &pc); CHKERRV(ierr);


  // check, if GAMG preconditioner is selected
  PetscBool useGAMGPreconditioner;
  PetscBool useHYPREPreconditioner;
  ierr = PetscObjectTypeCompare((PetscObject)pc, PCGAMG, &useGAMGPreconditioner); CHKERRV(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)pc, PCHYPRE, &useHYPREPreconditioner); CHKERRV(ierr);

  // only set coordinates ith the GAMG preconditioner is selected
  if ((useGAMGPreconditioner && FunctionSpaceType::Mesh::dim() < 3) || useHYPREPreconditioner)
  {
    // set the local node positions for the preconditioner
    int nDofsPerNode = data_.functionSpace()->nDofsPerNode();
    int nNodesLocal = data_.functionSpace()->nNodesLocalWithoutGhosts();

    std::vector<double> nodePositionCoordinatesForPreconditioner;
    nodePositionCoordinatesForPreconditioner.reserve(3*nNodesLocal);

    // loop over nodes and add their node positions
    for (dof_no_t dofNoLocal = 0; dofNoLocal < nNodesLocal*nDofsPerNode; dofNoLocal++)
    {
      Vec3 nodePosition = data_.functionSpace()->getGeometry(dofNoLocal);

      // add the coordinates
      for (int i = 0; i < 3; i++)
        nodePositionCoordinatesForPreconditioner.push_back(nodePosition[i]);
    }


    LOG(DEBUG) << "set coordinates to preconditioner, " << nodePositionCoordinatesForPreconditioner.size() << " node coordinates";
    ierr = PCSetCoordinates(pc, 3, nNodesLocal, nodePositionCoordinatesForPreconditioner.data()); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents,typename Term>
std::shared_ptr<typename FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::SlotConnectorDataType>
FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Term>::
getSlotConnectorData()
{
  return data_.getSlotConnectorData();
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents>
void FiniteElementMethodInitializeData<FunctionSpaceType,QuadratureType,nComponents,Equation::Dynamic::DirectionalDiffusion>::
initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction,
           std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> spatiallyVaryingPrefactor,
           bool useAdditionalDiffusionTensor)
{
  LOG(DEBUG) << "FiniteElementMethodInitializeData::initialize";

  // initialize the DiffusionTensorDirectional object
  this->data_.initialize(direction, spatiallyVaryingPrefactor, useAdditionalDiffusionTensor);

  // call normal initialize, this does not initialize the data object again
  FiniteElementMethodBase<FunctionSpaceType,QuadratureType,nComponents,Equation::Dynamic::DirectionalDiffusion>::initialize();
}

} // namespace

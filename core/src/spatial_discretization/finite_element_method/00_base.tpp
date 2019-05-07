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
#include "mesh/mesh_manager.h"
#include "solver/solver_manager.h"
#include "solver/linear.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

namespace SpatialDiscretization
{

template<typename FunctionSpaceType,typename QuadratureType,typename Term>
FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
FiniteElementMethodBase(DihuContext context, std::shared_ptr<FunctionSpaceType> functionSpace) :
  context_(context["FiniteElementMethod"]), data_(context["FiniteElementMethod"]), specificSettings_(context_.getPythonConfig()), initialized_(false)
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

template<typename FunctionSpaceType,typename QuadratureType,typename Term>
std::shared_ptr<FunctionSpaceType> FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
functionSpace()
{
  return data_.functionSpace();
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term>
Data::FiniteElements<FunctionSpaceType,Term> &FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
data()
{
  return data_;
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
setRankSubset(Partition::RankSubset rankSubset)
{
  data_.setRankSubset(rankSubset);
}
 
template<typename FunctionSpaceType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
initialize()
{
  // do not initialize if it was called already
  if (initialized_)
    return;

  data_.initialize();

  // assemble stiffness matrix
  Control::PerformanceMeasurement::start("durationSetStiffnessMatrix");
  setStiffnessMatrix();
  Control::PerformanceMeasurement::stop("durationSetStiffnessMatrix");

  // set the rhs
  Control::PerformanceMeasurement::start("durationSetRightHandSide");
  setRightHandSide();
  Control::PerformanceMeasurement::stop("durationSetRightHandSide");

  // apply boundary conditions
  Control::PerformanceMeasurement::start("durationAssembleBoundaryConditions");
  this->applyBoundaryConditions();
  Control::PerformanceMeasurement::stop("durationAssembleBoundaryConditions");

  initialized_ = true;
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
reset()
{
  data_.reset();
  initialized_ = false;
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
run()
{
  initialize();
  solve();
  data_.print();

  outputWriterManager_.writeOutput(data_);
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
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

  // non-zero initial values
#if 0  
  PetscScalar scalar = 0.5;
  ierr = VecSet(data_.solution()->values(), scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRV(ierr);
#endif


// get geometry field values of old mesh
  std::vector<Vec3> geometryFieldValues;
  std::vector<double> geometryFieldValuesConsecutive(geometryFieldValues.size()*3);

  this->data_.functionSpace()->geometryField().getValuesWithoutGhosts(geometryFieldValues);
  for (int i = 0; i < geometryFieldValues.size(); i++)
  {
    geometryFieldValuesConsecutive[3*i + 0] = geometryFieldValues[i][0];
    geometryFieldValuesConsecutive[3*i + 1] = geometryFieldValues[i][1];
    geometryFieldValuesConsecutive[3*i + 2] = geometryFieldValues[i][2];
  }
  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC(*ksp, &pc); CHKERRV(ierr);

  //ierr = PCSetCoordinates(pc, PetscInt dim, Pe tscInt nloc, geometryFieldValuesConsecutive.data());  CHKERRV(ierr);

  LOG(DEBUG) << "finiteElementMethod has geometryFieldValues: " << geometryFieldValues;

  LOG(DEBUG) << "solve...";

  // solve the system
  linearSolver->solve(data_.rightHandSide()->valuesGlobal(), data_.solution()->valuesGlobal(), "Solution obtained");
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
solveMG()
{
  
	if(kspMG_ != nullptr)
	{
		//if ksp object was already set up
		PetscErrorCode ierr;
    ierr = KSPSolve(*kspMG_, data_.rightHandSide()->valuesGlobal(),data_.solution()->valuesGlobal()); CHKERRV(ierr);
		int numberOfIterations = 0;
		PetscReal residualNorm = 0.0;
		ierr = KSPGetIterationNumber(*kspMG_, &numberOfIterations); CHKERRV(ierr);
		ierr = KSPGetResidualNorm(*kspMG_, &residualNorm); CHKERRV(ierr);

		KSPConvergedReason convergedReason;
		ierr = KSPGetConvergedReason(*kspMG_, &convergedReason); CHKERRV(ierr);

		LOG(INFO) << "Solution obtained in " << numberOfIterations << " iterations, residual norm " << residualNorm
		<< ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);
		} else {

		// solve linear system k*d=f for d
		LOG(TRACE) << "FiniteElementMethod::solveMG";

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

		// get linearMG solver context from solver manager
		std::shared_ptr<Solver::Linear> linearSolver = this->context_.solverManager()->template solver<Solver::Linear>(
		this->specificSettings_, this->data_.functionSpace()->meshPartition()->mpiCommunicator());
		std::shared_ptr<KSP> kspMG_ = linearSolver->ksp();

		assert(kspMG_ != nullptr);

		// set matrix used for linear system and preconditioner to ksp context
		PetscErrorCode ierr;
		ierr = KSPSetOperators(*kspMG_, stiffnessMatrix->valuesGlobal(), stiffnessMatrix->valuesGlobal()); CHKERRV(ierr);
		ierr =   KSPSetTolerances(*kspMG_,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1);

		// non-zero initial values
		#if 0  
		PetscScalar scalar = 0.5;
		ierr = VecSet(data_.solution()->values(), scalar); CHKERRV(ierr);
		ierr = KSPSetInitialGuessNonzero(*kspMG_, PETSC_TRUE); CHKERRV(ierr);

		#endif

		// solve the system
		ierr = KSPSolve(*kspMG_, data_.rightHandSide()->valuesGlobal(), data_.solution()->valuesGlobal()); CHKERRV(ierr);

		int numberOfIterations = 0;
        PetscReal residualNorm = 0.0;
		ierr = KSPGetIterationNumber(*kspMG_, &numberOfIterations); CHKERRV(ierr);
		ierr = KSPGetResidualNorm(*kspMG_, &residualNorm); CHKERRV(ierr);

		KSPConvergedReason convergedReason;
		ierr = KSPGetConvergedReason(*kspMG_, &convergedReason); CHKERRV(ierr);

		LOG(INFO) << "Solution obtained in " << numberOfIterations << " iterations, residual norm " << residualNorm
		<< ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);

		// check if solution is correct
		#if 0
		{
		// get rhs and solution from PETSc
		int vectorSize = 0;
		VecGetSize(data_.solution()->values(), &vectorSize);

		std::vector<int> indices(vectorSize);
		std::iota(indices.begin(), indices.end(), 0);
		std::vector<double> solution(vectorSize);
		std::vector<double> rhs(vectorSize);

		VecGetValues(data_.solution()->values(), vectorSize, indices.data(), solution.data());
		VecGetValues(data_.rightHandSide()->values(), vectorSize, indices.data(), rhs.data());

		// get stiffness matrix
		int nRows, nColumns;
		MatGetSize(data_.stiffnessMatrix(), &nRows, &nColumns);
		std::vector<int> rowIndices(nRows);
		std::iota(rowIndices.begin(), rowIndices.end(), 0);
		std::vector<int> columnIndices(nColumns);
		std::iota(columnIndices.begin(), columnIndices.end(), 0);
		std::vector<double> matrixValues(nRows*nColumns);

		std::vector<long int> nEntries = {nRows, nColumns};

		MatGetValues(data_.stiffnessMatrix().values(), nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());

		std::vector<double> f(vectorSize);

		// compute f = matrix * solution

		for (int i=0; i<vectorSize; i++)
		{
		  f[i] = 0.0;
		  for (int j=0; j<vectorSize; j++)
		  {
			f[i] += matrixValues[i*nColumns + j] * solution[j];
		  }
		}

		// compute residual norm
		double res = 0.0;
		for (int i=0; i<vectorSize; i++)
		{
		  res += (f[i] - rhs[i]) * (f[i] - rhs[i]);
		  LOG(DEBUG) << i << ". solution=" << solution[i]<< ", f=" <<f[i]<< ", rhs=" <<rhs[i]<< ", squared error: " <<(f[i] - rhs[i]) * (f[i] - rhs[i]);
		}

		LOG(DEBUG) << "res=" << res;
		}
		#endif 
	}
	
	outputWriterManager_.writeOutput(data_);

  VLOG(1) << "solution: " << *data_.solution();
}

template<typename FunctionSpaceType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
getResidual(Vec *residual)
{
    
    KSPBuildResidual(*kspMG_,NULL,NULL,residual);
}

template<typename FunctionSpaceType,typename QuadratureType>
void FiniteElementMethodInitializeData<FunctionSpaceType,QuadratureType,Equation::Dynamic::DirectionalDiffusion>::
initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction,
           std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> spatiallyVaryingPrefactor,
           bool useAdditionalDiffusionTensor)
{
  LOG(DEBUG) << "FiniteElementMethodInitializeData::initialize";

  // initialize the DiffusionTensorDirectional object
  this->data_.initialize(direction, spatiallyVaryingPrefactor, useAdditionalDiffusionTensor);

  // call normal initialize, this does not initialize the data object again
  FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Equation::Dynamic::DirectionalDiffusion>::initialize();
}

} // namespace

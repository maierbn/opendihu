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
#include "partition/partitioned_petsc_vec.h"
#include "partition/partitioned_petsc_mat.h"

namespace SpatialDiscretization
{

template<typename FunctionSpaceType,typename QuadratureType,typename Term>
FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Term>::
FiniteElementMethodBase(DihuContext context) :
  context_(context["FiniteElementMethod"]), data_(context["FiniteElementMethod"]), initialized_(false)
{
  specificSettings_ = context_.getPythonConfig();
  outputWriterManager_.initialize(context_, specificSettings_);

  // Create mesh or retrieve existing mesh from meshManager. This already creates meshPartition in functionSpace.initialize(), see function_space/03_function_space_partition_structured.tpp
  std::shared_ptr<Mesh::Mesh> mesh = context_.meshManager()->functionSpace<FunctionSpaceType>(specificSettings_);
  
  // store mesh in data
  data_.setFunctionSpace(std::static_pointer_cast<FunctionSpaceType>(mesh));
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
  setStiffnessMatrix();
  setRightHandSide();
  this->applyBoundaryConditions();

  initialized_ = true;
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

  // non-zero initial values
#if 0  
  PetscScalar scalar = 0.5;
  ierr = VecSet(data_.solution()->values(), scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRV(ierr);
#endif

  // solve the system
  ierr = KSPSolve(*ksp, data_.rightHandSide()->valuesGlobal(), data_.solution()->valuesGlobal()); CHKERRV(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp, &residualNorm); CHKERRV(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(*ksp, &convergedReason); CHKERRV(ierr);

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

    for(int i=0; i<vectorSize; i++)
    {
      f[i] = 0.0;
      for(int j=0; j<vectorSize; j++)
      {
        f[i] += matrixValues[i*nColumns + j] * solution[j];
      }
    }

    // compute residual norm
    double res = 0.0;
    for(int i=0; i<vectorSize; i++)
    {
      res += (f[i] - rhs[i]) * (f[i] - rhs[i]);
      LOG(DEBUG) << i << ". solution=" << solution[i]<< ", f=" <<f[i]<< ", rhs=" <<rhs[i]<< ", squared error: " <<(f[i] - rhs[i]) * (f[i] - rhs[i]);
    }

    LOG(DEBUG) << "res=" << res;
  }
#endif  
}

template<typename FunctionSpaceType,typename QuadratureType>
void FiniteElementMethodInitializeData<FunctionSpaceType,QuadratureType,Equation::Dynamic::DirectionalDiffusion>::
initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction, int multidomainNCompartments)
{
  LOG(DEBUG) << "FiniteElementMethodInitializeData::initialize";

  // initialize the DiffusionTensorFieldVariable object
  this->data_.initialize(direction, multidomainNCompartments);

  // call normal initialize, this does not initialize the data object again
  FiniteElementMethodBase<FunctionSpaceType,QuadratureType,Equation::Dynamic::DirectionalDiffusion>::initialize();
}

};

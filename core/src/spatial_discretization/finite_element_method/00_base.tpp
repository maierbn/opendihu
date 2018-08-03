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

template<typename BasisOnMeshType,typename QuadratureType,typename Term>
FiniteElementMethodBase<BasisOnMeshType,QuadratureType,Term>::
FiniteElementMethodBase(DihuContext context) :
  context_(context["FiniteElementMethod"]), data_(context["FiniteElementMethod"])
{
  specificSettings_ = context_.getPythonConfig();
  outputWriterManager_.initialize(specificSettings_);

  if (VLOG_IS_ON(2))
  {
    VLOG(2) << "FiniteElementMethodBase::FiniteElementMethodBase querying meshManager for mesh, specificSettings_:";
    PythonUtility::printDict(specificSettings_);
  }

  // Create mesh or retrieve existing mesh from meshManager. This does not yet create meshPartition, it is done later in data_.initialize().
  std::shared_ptr<Mesh::Mesh> mesh = context_.meshManager()->mesh<BasisOnMeshType>(specificSettings_);
  
  // store mesh in data
  data_.setMesh(std::static_pointer_cast<BasisOnMeshType>(mesh));
}

template<typename BasisOnMeshType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<BasisOnMeshType,QuadratureType,Term>::
applyBoundaryConditions()
{
  // PETSc Mat object for stiffness matrix needs to be assembled for this

  LOG(TRACE)<<"applyBoundaryConditions";

  dof_no_t nLocalUnknowns = this->data_.nLocalUnknowns();
  node_no_t nNodes = this->data_.mesh()->nLocalNodes();

  FieldVariable::FieldVariable<BasisOnMeshType,1> &rightHandSide = data_.rightHandSide();
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix = data_.stiffnessMatrix();
  
  // add Dirichlet boundary conditions
  // Boundary conditions are specified for dof numbers, not nodes, such that for Hermite it is possible to prescribe derivatives.
  // However the ordering of the dofs is not known in the config for unstructured meshes. Therefore the ordering is special.
  // For every node there are as many values as dofs, in contiguous order.
  // Example for 2D Hermite, unstructured grid, 2x2 elements:
  //
  // node numbering:
  //  6_7_8
  // 3|_4_|5
  // 0|_1_|2
  //
  // dof numbering:
  //  6_7_8
  // 2|_3_|5
  // 0|_1_|4
  //
  // To specify du/dn = 0 an the left boundary you would set:
  // bc[0*2+1] = 0, bc[3*2+1] = 0, bc[6*2+1] = 0
  //
  // To specifiy u=0 on the bottom, you would set:
  // bc[0] = 0, bc[2] = 0, bc[4] = 0
  

  // get the first dirichlet boundary condition from the list
  std::pair<node_no_t, double> boundaryCondition
    = PythonUtility::getOptionDictBegin<node_no_t, double>(specificSettings_, "DirichletBoundaryCondition");

  // loop over Dirichlet boundary conditions
  for (; !PythonUtility::getOptionDictEnd(specificSettings_, "DirichletBoundaryCondition");
       PythonUtility::getOptionDictNext<node_no_t, double>(specificSettings_, "DirichletBoundaryCondition", boundaryCondition))
  {
    dof_no_t boundaryConditionIndex = boundaryCondition.first;
    double boundaryConditionValue = boundaryCondition.second;

    // omit negative indices
    if (boundaryConditionIndex < 0)
      continue;

    // translate BC index to nodeNo and dofIndex
    node_no_t boundaryConditionNodeNo = boundaryConditionIndex / BasisOnMeshType::nDofsPerNode();
    int boundaryConditionNodalDofIndex = boundaryConditionIndex - boundaryConditionNodeNo * BasisOnMeshType::nDofsPerNode();
    
    if (boundaryConditionIndex > nLocalUnknowns)
    {
      LOG(WARNING) << "Boundary condition specified for index " << boundaryConditionIndex
        << " (on local node " << boundaryConditionNodeNo << ", index " << boundaryConditionNodalDofIndex << ")"
        << ", but scenario has only " << nLocalUnknowns << " local unknowns, " << nNodes << " local nodes";
       continue;
    }
    
    //TODO handle if boundary conditions are given as global indices
    
    dof_no_t boundaryConditionDofNo = this->data_.mesh()->getNodeDofNo(boundaryConditionNodeNo, boundaryConditionNodalDofIndex);

    // set rhs entry to prescribed value
    rightHandSide.setValue(boundaryConditionDofNo, boundaryConditionValue, INSERT_VALUES);

    VLOG(1) << "  BC node " << boundaryConditionNodeNo << " index " << boundaryConditionIndex 
      << ", dof " << boundaryConditionDofNo << ", value " << boundaryConditionValue;

    // get the column number boundaryConditionDofNo of the stiffness matrix. It is needed for updating the rhs.
    std::vector<int> rowIndices((int)nLocalUnknowns);
    std::iota (rowIndices.begin(), rowIndices.end(), 0);    // fill with increasing numbers: 0,1,2,...
    std::vector<int> columnIndices = {(int)boundaryConditionDofNo};

    std::vector<double> coefficients(nLocalUnknowns);

    stiffnessMatrix->getValuesGlobalIndexing(nLocalUnknowns, rowIndices.data(), 1, columnIndices.data(), coefficients.data());

    // set values of row and column of the DOF to zero and diagonal entry to 1
    int matrixIndex = (int)boundaryConditionDofNo;
    stiffnessMatrix->zeroRowsColumns(1, &matrixIndex, 1.0);
    
    // update rhs
    for (node_no_t rowNo = 0; rowNo < nLocalUnknowns; rowNo++)
    {
      if (rowNo == boundaryConditionDofNo)
       continue;

      // update rhs value to be f_new = f_old - m_{ij}*u_{i} where i is the index of the prescribed node,
      // m_{ij} is entry of stiffness matrix and u_{i} is the prescribed value
      double rhsSummand = -coefficients[rowNo] * boundaryConditionValue;
      rightHandSide.setValue(rowNo, rhsSummand, ADD_VALUES);

      LOG_IF(false,DEBUG) << "  in row " << rowNo << " add " << rhsSummand << " to rhs, coefficient: " << coefficients[rowNo];
    }
  }
}

template<typename BasisOnMeshType,typename QuadratureType,typename Term>
std::shared_ptr<Mesh::Mesh> FiniteElementMethodBase<BasisOnMeshType,QuadratureType,Term>::
mesh()
{
  return data_.mesh();
}

template<typename BasisOnMeshType,typename QuadratureType,typename Term>
Data::FiniteElements<BasisOnMeshType,Term> &FiniteElementMethodBase<BasisOnMeshType,QuadratureType,Term>::
data()
{
  return data_;
}

template<typename BasisOnMeshType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<BasisOnMeshType,QuadratureType,Term>::
setRankSubset(Partition::RankSubset rankSubset)
{
  data_.setRankSubset(rankSubset);
}
 
template<typename BasisOnMeshType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<BasisOnMeshType,QuadratureType,Term>::
initialize()
{
  data_.initialize();
  setStiffnessMatrix();
  setRightHandSide();
  data_.finalAssembly();
  applyBoundaryConditions();
}

template<typename BasisOnMeshType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<BasisOnMeshType,QuadratureType,Term>::
run()
{
  initialize();
  solve();
  data_.print();

  outputWriterManager_.writeOutput(data_);
}

template<typename BasisOnMeshType,typename QuadratureType,typename Term>
void FiniteElementMethodBase<BasisOnMeshType,QuadratureType,Term>::
solve()
{
  // solve linear system k*d=f for d
  LOG(TRACE) << "FiniteElementMethod::solve";

  // if equation was set to none, do not solve the problem (this is for unit tests that don't test for solution)
  if (std::is_same<Term,Equation::None>::value)
   return;

  // get stiffness matrix
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix = data_.stiffnessMatrix();

  // get linear solver context from solver manager
  std::shared_ptr<Solver::Linear> linearSolver = this->context_.solverManager()->template solver<Solver::Linear>(this->specificSettings_);
  std::shared_ptr<KSP> ksp = linearSolver->ksp();
  assert(ksp != nullptr);

  // set matrix used for linear system and preconditioner to ksp context
  PetscErrorCode ierr;
  ierr = KSPSetOperators (*ksp, stiffnessMatrix->values(), stiffnessMatrix->values()); CHKERRV(ierr);

  // non-zero initial values
#if 0  
  PetscScalar scalar = 0.5;
  ierr = VecSet(data_.solution().values(), scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRV(ierr);
#endif

  // solve the system
  ierr = KSPSolve(*ksp, data_.rightHandSide().values(), data_.solution().values()); CHKERRV(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp, &residualNorm); CHKERRV(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(*ksp, &convergedReason); CHKERRV(ierr);

  LOG(INFO) << "Solution done in " << numberOfIterations << " iterations, residual norm " << residualNorm
    << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);

  // check if solution is correct
#if 0
  {
    // get rhs and solution from PETSc
    int vectorSize = 0;
    VecGetSize(data_.solution().values(), &vectorSize);

    std::vector<int> indices(vectorSize);
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<double> solution(vectorSize);
    std::vector<double> rhs(vectorSize);

    VecGetValues(data_.solution().values(), vectorSize, indices.data(), solution.data());
    VecGetValues(data_.rightHandSide().values(), vectorSize, indices.data(), rhs.data());

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
      LOG(DEBUG) << i << ". solution="<<solution[i]<<", f="<<f[i]<<", rhs="<<rhs[i]<<", squared error: "<<(f[i] - rhs[i]) * (f[i] - rhs[i]);
    }

    LOG(DEBUG) << "res=" << res;
  }
#endif  
}

};
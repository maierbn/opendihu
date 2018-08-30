#include "spatial_discretization/finite_element_method/02_boundary_conditions.h"

#include <Python.h>
#include <memory>
#include <vector>
#include <petscsys.h>

#include "quadrature/tensor_product.h"

namespace SpatialDiscretization
{

template<typename FunctionSpaceType,typename Term>
void BoundaryConditions<FunctionSpaceType,Quadrature::None,Term,Term>::
applyBoundaryConditions()
{
  applyBoundaryConditionsStrongForm();
}

template<typename FunctionSpaceType,typename Term>
void BoundaryConditions<FunctionSpaceType,Quadrature::None,Term,Term>::
applyBoundaryConditionsStrongForm()
{
  // This sets rows and columns in stiffness matrix to 0 and diagonal to 1 for BC dofs.
  // This only works in serial, because getValuesGlobalIndexing does not work for distributed matrices
  // PETSc Mat object for stiffness matrix needs to be assembled for this.

  LOG(TRACE) << "applyBoundaryConditionsStrongForm";

  dof_no_t nDofsLocal = this->data_.functionSpace()->nDofsLocalWithoutGhosts();
  node_no_t nNodes = this->data_.functionSpace()->nNodesLocalWithoutGhosts();

  FieldVariable::FieldVariable<FunctionSpaceType,1> &rightHandSide = this->data_.rightHandSide();
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix = this->data_.stiffnessMatrix();

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
  = PythonUtility::getOptionDictBegin<node_no_t, double>(this->specificSettings_, "DirichletBoundaryCondition");

  // loop over Dirichlet boundary conditions
  for (; !PythonUtility::getOptionDictEnd(this->specificSettings_, "DirichletBoundaryCondition");
  PythonUtility::getOptionDictNext<node_no_t, double>(this->specificSettings_, "DirichletBoundaryCondition", boundaryCondition))
  {
    dof_no_t boundaryConditionIndex = boundaryCondition.first;
    double boundaryConditionValue = boundaryCondition.second;

    // omit negative indices
    if (boundaryConditionIndex < 0)
      continue;

    // translate BC index to nodeNo and dofIndex
    node_no_t boundaryConditionNodeNo = boundaryConditionIndex / FunctionSpaceType::nDofsPerNode();
    int boundaryConditionNodalDofIndex = boundaryConditionIndex - boundaryConditionNodeNo * FunctionSpaceType::nDofsPerNode();

    if (boundaryConditionIndex > nDofsLocal)
    {
      LOG(WARNING) << "Boundary condition specified for index " << boundaryConditionIndex
      << " (on local node " << boundaryConditionNodeNo << ", index " << boundaryConditionNodalDofIndex << ")"
      << ", but scenario has only " << nDofsLocal << " local unknowns, " << nNodes << " local nodes";
      continue;
    }

    dof_no_t boundaryConditionDofNo = this->data_.functionSpace()->getNodeDofNo(boundaryConditionNodeNo, boundaryConditionNodalDofIndex);

    // set rhs entry to prescribed value
    rightHandSide.setValue(boundaryConditionDofNo, boundaryConditionValue, INSERT_VALUES);

    VLOG(1) << "  BC node " << boundaryConditionNodeNo << " index " << boundaryConditionIndex
    << ", dof " << boundaryConditionDofNo << ", value " << boundaryConditionValue;

    // get the column number boundaryConditionDofNo of the stiffness matrix. It is needed for updating the rhs.
    std::vector<int> rowIndices((int)nDofsLocal);
    std::iota (rowIndices.begin(), rowIndices.end(), 0);    // fill with increasing numbers: 0,1,2,...
    std::vector<int> columnIndices = {(int)boundaryConditionDofNo};

    std::vector<double> coefficients(nDofsLocal);

    stiffnessMatrix->getValuesGlobalIndexing(nDofsLocal, rowIndices.data(), 1, columnIndices.data(), coefficients.data());

    // set values of row and column of the DOF to zero and diagonal entry to 1
    int matrixIndex = (int)boundaryConditionDofNo;
    stiffnessMatrix->assembly(MAT_FINAL_ASSEMBLY);
    stiffnessMatrix->zeroRowsColumns(1, &matrixIndex, 1.0);

    // update rhs
    for (node_no_t rowNo = 0; rowNo < nDofsLocal; rowNo++)
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

};  // namespace

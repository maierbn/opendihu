#include "data_management/finite_elements.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"
#include "basis_on_mesh/basis_on_mesh.h"
#include "mesh/unstructured_deformable.h"
#include "basis_function/hermite.h"
#include "partition/partitioned_petsc_mat.h"

namespace Data
{

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
FiniteElements(DihuContext context) : Data<BasisOnMeshType>(context)
{
  LOG(TRACE) << "Data::FiniteElements constructor";
  //PythonUtility::printDict(this->context_.getPythonConfig());
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
~FiniteElements()
{
  // free PETSc objects
  //if (this->initialized_)
  //{
    //ierr = MatDestroy(&this->stiffnessMatrix_); CHKERRV(ierr);
  //}
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
initialize()
{
  Data<BasisOnMeshType>::initialize();

  // set up diffusion tensor if there is any
  DiffusionTensor<BasisOnMeshType::dim()>::initialize(this->context_.getPythonConfig());
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros)
{
  const int D = this->mesh_->dimension();
  const int nDofsPerNode = BasisOnMesh::BasisOnMeshBaseDim<1,typename BasisOnMeshType::BasisFunction>::nDofsPerNode();
  const int nDofsPerBasis = BasisOnMesh::BasisOnMeshBaseDim<1,typename BasisOnMeshType::BasisFunction>::nDofsPerElement();
  const int nOverlaps = (nDofsPerBasis*2 - 1) * nDofsPerNode;   // number of nodes of 2 neighbouring 1D elements (=number of ansatz functions in support of center ansatz function)

  // due to PETSc storage diagonalNonZeros and offdiagonalNonZeros should be both set to the maximum number of non-zero entries per row

  switch (D)
  {
  case 1:
    diagonalNonZeros = nOverlaps;
    offdiagonalNonZeros = diagonalNonZeros;
    break;
  case 2:
    diagonalNonZeros = pow(nOverlaps, 2) + 16;   // because of boundary conditions there can be more entries, which are all zero, but stored as non-zero
    offdiagonalNonZeros = diagonalNonZeros;
    break;
  case 3:
    diagonalNonZeros = pow(nOverlaps, 3);
    offdiagonalNonZeros = diagonalNonZeros;
    break;
  };
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
createPetscObjects()
{
  LOG(DEBUG)<<"FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::createPetscObjects";

  // get the partitioning from the mesh
  std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition = this->mesh_->meshPartition();
  
  // create field variables on local partition
  this->rhs_ = this->mesh_->template createFieldVariable<1>("rhs");
  this->solution_ = this->mesh_->template createFieldVariable<1>("solution");

  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);

  LOG(DEBUG) << "d="<<this->mesh_->dimension()
    <<", number of diagonal non-zeros: "<<diagonalNonZeros<<", number of off-diagonal non-zeros: "<<offdiagonalNonZeros;

  int nComponents = 1;
  this->stiffnessMatrix_ = std::make_shared<PartitionedPetscMat<BasisOnMeshType>>(meshPartition, nComponents, diagonalNonZeros, offdiagonalNonZeros, "stiffnessMatrix");
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
finalAssembly()
{
  this->stiffnessMatrix_->assembly(MAT_FINAL_ASSEMBLY);
  
  LOG(DEBUG) << "finalAssembly";
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
stiffnessMatrix()
{
  return this->stiffnessMatrix_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
FieldVariable::FieldVariable<BasisOnMeshType,1> &FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
rightHandSide()
{
  return *this->rhs_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
FieldVariable::FieldVariable<BasisOnMeshType,1> &FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
solution()
{
  return *this->solution_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
massMatrix()
{
  return this->massMatrix_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->stiffnessMatrix_;
  VLOG(4) << *this->rhs_;
  VLOG(4) << *this->solution_;
  
  MatInfo info;
  MatGetInfo(this->stiffnessMatrix_->valuesGlobal(), MAT_LOCAL, &info);

  VLOG(4)<<"stiffnessMatrix info: "<<std::endl
    <<"block_size: "<<info.block_size<<std::endl
    <<"number of nonzeros: allocated: "<<info.nz_allocated<<", used: "<<info.nz_used<<", unneeded: "<<info.nz_unneeded<<std::endl
    <<"memory allocated: "<<info.memory<<std::endl
    <<"number of matrix assemblies called: "<<info.assemblies<<std::endl
    <<"number of mallocs during MatSetValues(): "<<info.mallocs<<std::endl
    <<"fill ratio for LU/ILU: given: "<<info.fill_ratio_given<<", needed: "<<info.fill_ratio_needed<<std::endl
    <<"number of mallocs during factorization: "<<info.factor_mallocs<<std::endl;

  VLOG(4)<<"======================";
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
bool FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
massMatrixInitialized()
{
  return this->massMatrixInitialized_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
initializeMassMatrix()
{
  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);

  std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> partition = this->mesh_->meshPartition();
  const int nComponents = 1;
  this->massMatrix_ = std::make_shared<PartitionedPetscMat<BasisOnMeshType>>(partition, nComponents, diagonalNonZeros, offdiagonalNonZeros, "massMatrix");
  
  this->massMatrixInitialized_ = true;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
typename FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::OutputFieldVariables FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
getOutputFieldVariables()
{
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<BasisOnMeshType,3>>(this->mesh_->geometryField());
  /*
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> generalField;

  generalField = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,3>>(this->mesh_->fieldVariable("general"));
  if (!generalField)
   generalField = geometryField;
  */
  return OutputFieldVariables(
    geometryField,
    solution_,
    rhs_
  );
}


} // namespace Data

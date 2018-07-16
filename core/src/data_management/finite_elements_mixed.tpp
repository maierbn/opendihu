#include "data_management/finite_elements_mixed.h"

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

namespace Data
{

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
FieldVariable::FieldVariable<LowOrderBasisOnMeshType,1> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
pressure()
{
  return *this->pressure_;
}


template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
typename FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::OutputFieldVariables FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
getOutputFieldVariables()
{
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,3>> actualGeometryField
    = std::make_shared<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,3>>(this->mesh_->geometryField());

  return OutputFieldVariables(
    this->geometryReference_,
    actualGeometryField,
    this->displacements_,
    this->pressure_,
    this->residual_,
    this->externalVirtualWork_
  );
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
setMesh(std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> mixedMesh)
{
  mixedMesh_ = mixedMesh;

  // store high order mesh as mesh_
  this->mesh_ = mixedMesh_->highOrderBasisOnMesh();
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
mixedMesh()
{
  return mixedMesh_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
initialize()
{
  // call initialize of parent class
  FiniteElementsSolidMechanics<HighOrderBasisOnMeshType,Term>::initialize();

  LOG(DEBUG) << "mesh has geometry field: " << this->mesh_->hasGeometryField();
  initializeFieldVariables();
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
initializeFieldVariables()
{
  std::vector<std::string> unnamedSingleComponent({"0"});
  this->pressure_ = std::make_shared<FieldVariable::FieldVariable<LowOrderBasisOnMeshType,1>>(this->mixedMesh_->lowOrderBasisOnMesh(), "pressure", unnamedSingleComponent);
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
const dof_no_t FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
getTangentStiffnessMatrixNRows()
{
  const int D = HighOrderBasisOnMeshType::dim();
  return this->mixedMesh_->highOrderBasisOnMesh()->nLocalDofs() * D + this->mixedMesh_->lowOrderBasisOnMesh()->nLocalDofs();
}

} // namespace Data

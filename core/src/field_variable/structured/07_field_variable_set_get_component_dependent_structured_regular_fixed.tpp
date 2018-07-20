#include "field_variable/structured/07_field_variable_set_get_component_dependent_structured_regular_fixed.h"

#include <sstream>
#include "easylogging++.h"
#include "utility/string_utility.h"

namespace FieldVariable
{

//! get a single value from global dof no. for all components
template<int D,typename BasisFunctionType,int nComponents>
std::array<double,nComponents> FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>::
getValue(node_no_t dofGlobalNo)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    std::array<PetscInt,nComponents> indices;
    std::array<double,nComponents> result;

    // get number of dofs
    assert(this->mesh_);
    const dof_no_t nDofs = this->mesh_->nLocalDofs();

    // prepare lookup indices for PETSc vector values_
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      indices[componentNo] = componentNo*nDofs + dofGlobalNo;
    }

    VecGetValues(this->values_, nComponents, indices.data(), result.data());
    return result;
  }

  // if this is a geometry field compute the information
  const node_no_t nNodesInXDirection = this->mesh_->nNodes(0);
  const node_no_t nNodesInYDirection = this->mesh_->nNodes(1);
  const int nDofsPerNode = BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::nDofsPerNode();

  int nodeNo = int(dofGlobalNo / nDofsPerNode);
  int nodeLocalDofIndex = int(dofGlobalNo % nDofsPerNode);

  std::array<double,nComponents> value;
  if (nodeLocalDofIndex > 0)   // if this is a derivative of Hermite, set to 0
  {
    value[0] = 0;
    value[1] = 0;
    value[2] = 0;
  }
  else
  {
    // x direction
    value[0] = (nodeNo % nNodesInXDirection) * this->meshWidth_;

    // y direction
    value[1] = (int(nodeNo / nNodesInXDirection) % nNodesInYDirection) * this->meshWidth_;

    // z direction
    value[2] = int(nodeNo / (nNodesInXDirection*nNodesInYDirection)) * this->meshWidth_;
  }
  return value;
}

};

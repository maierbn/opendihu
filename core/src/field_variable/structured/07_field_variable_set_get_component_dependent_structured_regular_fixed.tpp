#include "field_variable/structured/07_field_variable_set_get_component_dependent_structured_regular_fixed.h"

#include <sstream>
#include "easylogging++.h"
#include "utility/string_utility.h"

namespace FieldVariable
{

//! get a single value from local dof no. for all components
template<int D,typename BasisFunctionType,int nComponents>
std::array<double,nComponents> FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>::
getValue(node_no_t dofLocalNo) const
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    std::array<double,nComponents> result;

    // prepare lookup indices for PETSc vector values_
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      this->values_->getValues(componentNo, 1, &dofLocalNo, result.data()+componentNo);
    }

    return result;
  }

  assert(dofLocalNo < this->functionSpace_->nDofsWithoutGhosts());
   
  // for geometry field compute the entries, this does not work for ghost dofs
  const node_no_t nLocalNodesInXDirection = this->functionSpace_->nNodesLocalWithoutGhosts(0);
  const node_no_t nLocalNodesInYDirection = this->functionSpace_->nNodesLocalWithoutGhosts(1);
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::nDofsPerNode();
  
  int nodeLocalNo = int(dofLocalNo / nDofsPerNode);
  int nodeLocalDofIndex = int(dofLocalNo % nDofsPerNode);

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
    value[0] = this->functionSpace_->meshPartition()->beginNodeGlobal(0) * this->functionSpace_->meshWidth()
      + (nodeLocalNo % nLocalNodesInXDirection) * this->functionSpace_->meshWidth();

    // y direction
    value[1] = this->functionSpace_->meshPartition()->beginNodeGlobal(1) * this->functionSpace_->meshWidth()
      + (int(nodeLocalNo / nLocalNodesInXDirection) % nLocalNodesInYDirection) * this->functionSpace_->meshWidth();

    // z direction
    value[2] = this->functionSpace_->meshPartition()->beginNodeGlobal(2) * this->functionSpace_->meshWidth()
      + int(nodeLocalNo / (nLocalNodesInXDirection*nLocalNodesInYDirection)) * this->functionSpace_->meshWidth();
  }
  
  return value;
}

};

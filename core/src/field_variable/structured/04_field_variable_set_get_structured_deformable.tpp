#include "field_variable/structured/04_field_variable_set_get_structured_deformable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <cassert>

namespace FieldVariable
{
 

//! copy the values from another field variable of the same type
template<int D,typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
setValues(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &rhs)
{
  VecCopy(rhs.values_, this->values_);
}

};
#include "field_variable/structured/05_field_variable_data_structured_regular_fixed.h"

#include <sstream>
#include "easylogging++.h"
#include "utility/string_utility.h"

namespace FieldVariable
{

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>::
setMeshWidth(double meshWidth)
{
  meshWidth_ = meshWidth;
}

template<int D, typename BasisFunctionType, int nComponents>
double FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>::
meshWidth() const
{
  return this->meshWidth_;
}

};

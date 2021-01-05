#include "function_space/function_space.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "basis_function/lagrange.h"
#include "basis_function/hermite.h"
#include "mesh/structured_deformable.h"
#include "mesh/structured_regular_fixed.h"

namespace FunctionSpace
{

template<typename MeshType, typename BasisFunctionType>
std::string FunctionSpace<MeshType,BasisFunctionType>::
getDescription() const
{
  const int D = MeshType::dim();
  std::stringstream description;
  description << "\"" << this->meshName() << "\", " << D << "D ";

  if (std::is_same<MeshType, ::Mesh::StructuredDeformableOfDimension<D>>::value)
  {
    description << "structured deformable, ";
  }
  else if (std::is_same<MeshType, ::Mesh::StructuredRegularFixedOfDimension<D>>::value)
  {
    description << "regular fixed, ";
  }
  else if (std::is_same<MeshType, ::Mesh::UnstructuredDeformableOfDimension<D>>::value)
  {
    description << "unstructured, ";
  }

  if (std::is_same<BasisFunctionType, ::BasisFunction::LagrangeOfOrder<1>>::value)
  {
    description << "linear Lagrange basis";
  }
  else if (std::is_same<BasisFunctionType, ::BasisFunction::LagrangeOfOrder<2>>::value)
  {
    description << "quadratic Lagrange basis";
  }
  else if (std::is_same<BasisFunctionType, ::BasisFunction::Hermite>::value)
  {
    description << "cubic Hermite basis";
  }
  return description.str();
}

} // namespace

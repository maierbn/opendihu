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
std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> FunctionSpace<MeshType,BasisFunctionType>::
getElementDofNosLocal(element_no_t elementNo) const
{
  std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> dofNosLocal;
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    dofNosLocal[dofIndex] = this->getDofNo(elementNo, dofIndex);
  }
  return dofNosLocal;
}

//! vectorized version of getElementDofNosLocal
template<typename MeshType, typename BasisFunctionType>
std::array<Vc::int_v,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> FunctionSpace<MeshType,BasisFunctionType>::
getElementDofNosLocal(Vc::int_v elementNo) const
{
  std::array<Vc::int_v,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> dofNosLocal;
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    for (int vcComponent = 0; vcComponent < Vc::double_v::size(); vcComponent++)
    {
      if (elementNo[vcComponent] == -1)
        dofNosLocal[dofIndex][vcComponent] = -1;
      else
        dofNosLocal[dofIndex][vcComponent] = this->getDofNo(elementNo[vcComponent], dofIndex);
    }
  }
  return dofNosLocal;
}

template<typename MeshType, typename BasisFunctionType>
void FunctionSpace<MeshType,BasisFunctionType>::
getElementDofNosLocal(element_no_t elementNo, std::vector<dof_no_t> &dofNosLocal) const
{
  dofNosLocal.resize(FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement());
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    dofNosLocal[dofIndex] = this->getDofNo(elementNo, dofIndex);
  }
}

template<typename MeshType, typename BasisFunctionType>
void FunctionSpace<MeshType,BasisFunctionType>::
getElementDofNosLocalWithoutGhosts(element_no_t elementNo, std::vector<dof_no_t> &dofNosLocal) const
{
  dofNosLocal.reserve(FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement());
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    dof_no_t dofNoLocal = this->getDofNo(elementNo, dofIndex);
    if (dofNoLocal < this->meshPartition_->nDofsLocalWithoutGhosts())
    {
      dofNosLocal.push_back(dofNoLocal);
    }
  }
}

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

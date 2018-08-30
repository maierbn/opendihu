#include "function_space/function_space.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace FunctionSpace
{
template<typename MeshType,int D,int order>
std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> FunctionSpace<MeshType,BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>::
getElementDofNosLocal(element_no_t elementNo) const
{
  const int nDofsPerElement = FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement();
  std::array<dof_no_t,nDofsPerElement> dof;

  // fill array with increasing values, starting with elementNo*nDofsPerElement
  std::iota(dof.begin(), dof.end(), elementNo*nDofsPerElement);

  return dof;
}

template<typename MeshType,int D,int order>
void FunctionSpace<MeshType,BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>::
getElementDofNosLocal(element_no_t elementNo, std::vector<dof_no_t> &globalDofNos) const
{
  const int nDofsPerElement = FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement();
  globalDofNos.resize(nDofsPerElement);

  // fill array with increasing values, starting with elementNo*nDofsPerElement
  std::iota(globalDofNos.begin(), globalDofNos.end(), elementNo*nDofsPerElement);
}

};  // namespace
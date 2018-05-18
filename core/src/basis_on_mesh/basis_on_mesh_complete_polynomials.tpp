#include "basis_on_mesh/basis_on_mesh.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace BasisOnMesh
{
template<typename MeshType,int D,int order>
std::array<dof_no_t,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> BasisOnMesh<MeshType,BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>::
getElementDofNos(element_no_t elementNo) const
{
  const int nDofsPerElement = BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement();
  std::array<dof_no_t,nDofsPerElement> dof;

  // fill array with increasing values, starting with elementNo*nDofsPerElement
  std::iota(dof.begin(), dof.end(), elementNo*nDofsPerElement);

  return dof;
}

template<typename MeshType,int D,int order>
void BasisOnMesh<MeshType,BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>::
getElementDofNos(element_no_t elementNo, std::vector<dof_no_t> &globalDofNos) const
{
  const int nDofsPerElement = BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement();
  globalDofNos.resize(nDofsPerElement);

  // fill array with increasing values, starting with elementNo*nDofsPerElement
  std::iota(globalDofNos.begin(), globalDofNos.end(), elementNo*nDofsPerElement);
}

};  // namespace
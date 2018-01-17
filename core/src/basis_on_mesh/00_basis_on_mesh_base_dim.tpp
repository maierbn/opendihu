#include "basis_on_mesh/00_basis_on_mesh_base_dim.h"

#include "utility/math_utility.h"
#include "easylogging++.h"

#include <cmath>
#include <array>

namespace BasisOnMesh
{
  
template<int D,typename BasisFunctionType>
constexpr int BasisOnMeshBaseDim<D,BasisFunctionType>::
nDofsPerElement()
{
  return pow(BasisFunctionType::nDofsPerBasis(),D);
}

template<int D,typename BasisFunctionType>
constexpr int BasisOnMeshBaseDim<D,BasisFunctionType>::
nNodesPerElement()
{
  return pow(BasisFunctionType::nDofsPerBasis()/BasisFunctionType::nDofsPerNode(),D);
}

template<int D,typename BasisFunctionType>
constexpr int BasisOnMeshBaseDim<D,BasisFunctionType>::
nDofsPerNode()
{
  return nDofsPerElement() / nNodesPerElement();
}

template<int D,typename BasisFunctionType>
constexpr int BasisOnMeshBaseDim<D,BasisFunctionType>::
averageNNodesPerElement()
{
  return pow(BasisFunctionType::nDofsPerBasis()/BasisFunctionType::nDofsPerNode()-1,D);
}

template<int D,typename BasisFunctionType>
constexpr int BasisOnMeshBaseDim<D,BasisFunctionType>::
averageNDofsPerElement()
{
 // nNodesPerBasis = nDofsPerBasis / nDofsPerNode
 // averageNNodesPerBasis = nNodesPerBasis - 1
 // averageNDofsPerElement = averageNNodesPerBasis * nDofsPerNode
  return BasisOnMeshBaseDim<D,BasisFunctionType>::averageNNodesPerElement() * BasisFunctionType::nDofsPerNode();
}

};  // namespace
#include "basis_on_mesh/00_basis_on_mesh_base_dim.h"

#include "utility/math_utility.h"
#include "easylogging++.h"
#include "basis_function/hermite.h"
#include "basis_function/lagrange.h"

#include <cmath>
#include <array>

namespace BasisOnMesh
{
 
// lagrange and hermite polynomials  
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

// -------------------------
// complete polynomials 
template<int D,int order>
constexpr int BasisOnMeshBaseDim<D,BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>>::
nDofsPerElement()
{
  return BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>::nDofsPerBasis();
}
 
template<int D,int order>
constexpr int BasisOnMeshBaseDim<D,BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>>::
nNodesPerElement()
{
  return 0;
}

template<int D,int order>
constexpr int BasisOnMeshBaseDim<D,BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>>::
nDofsPerNode()
{
  return 0;
}

template<int D,int order>
constexpr int BasisOnMeshBaseDim<D,BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>>::
averageNDofsPerElement()
{
  return nDofsPerElement();
}

template<int D,int order>
constexpr int BasisOnMeshBaseDim<D,BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>>::
averageNNodesPerElement()
{
  return 0;
}

};  // namespace
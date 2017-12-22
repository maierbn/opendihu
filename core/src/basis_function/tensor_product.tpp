#include "basis_function/tensor_product.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace BasisFunction
{
  
template<int D,typename BasisFunctionType>
std::array<int,TensorProductBase<D,BasisFunctionType>::nDofsPerElement()> TensorProduct<D,BasisFunctionType>::
getElementDofs(element_idx_t elementNo, std::array<int,D> nElements)
{
  std::array<int,TensorProductBase<D,BasisFunctionType>::nDofsPerElement()> dof;
  for (int dofIndex = 0; dofIndex < TensorProductBase<D,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    dof[dofIndex] = TensorProductDim<D,BasisFunctionType>::getDofNo(elementNo, dofIndex, nElements);
  }
  
#if 0
  std::stringstream s;
  for (int i=0; i<D; i++)
  {
    s << nElements[i] << " ";
  }
  LOG(DEBUG) << "element " << elementNo << ", nElements: " << s.str() << ", nDofs: " << TensorProductBase<D,BasisFunctionType>::nDofsPerElement();
  s.str("");
  for (int dofIndex = 0; dofIndex < TensorProductBase<D,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    s << dof[dofIndex] << " ";
  }
    
  LOG(DEBUG) << "dofs: " << s.str();
#endif

  return dof;
}

template<int D,typename BasisFunctionType>
std::array<int,TensorProductBase<D,BasisFunctionType>::nNodesPerElement()> TensorProduct<D,BasisFunctionType>::
getElementNodes(element_idx_t elementNo, std::array<int,D> nElements)
{
  std::array<int,TensorProductBase<D,BasisFunctionType>::nNodesPerElement()> nodes;
  for (int nodeIndex = 0; nodeIndex < TensorProductBase<D,BasisFunctionType>::nNodesPerElement(); nodeIndex++)
  {
    nodes[nodeIndex] = TensorProductDim<D,BasisFunctionType>::getNodeNo(elementNo, nodeIndex, nElements);
  }
  return nodes;
}

template<int D,typename BasisFunctionType>
std::array<std::array<double,D>,TensorProductBase<D,BasisFunctionType>::nDofsPerElement()> TensorProduct<D,BasisFunctionType>::
getGradPhi(std::array<double,D> xi)
{
  std::array<std::array<double,D>,TensorProductBase<D,BasisFunctionType>::nDofsPerElement()> gradPhi;
  for (int dofIndex = 0; dofIndex < TensorProductBase<D,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    gradPhi[dofIndex] = TensorProductDim<D,BasisFunctionType>::gradPhi(dofIndex, xi);
  }
  return gradPhi;
}

};  // namespace
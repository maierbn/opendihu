#include "function_space/09_function_space_find_position.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "mesh/face_t.h"

namespace FunctionSpace
{

// unstructured
template<int D,typename BasisFunctionType>
bool FunctionSpaceFindPosition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType,Mesh::UnstructuredDeformableOfDimension<D>>::
findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,D> &xi, bool startSearchInCurrentElement, double &residual, bool &searchedAllElements, double xiTolerance)
{
  const element_no_t nElements = this->nElementsLocal();

  // set starting no to 0 if it was not given and is thus arbitrarily initialized
  if (elementNo < 0 || elementNo >= nElements)
    elementNo = 0;

  // check if point is already in current element
  if (this->pointIsInElement(point, elementNo, xi, residual, xiTolerance))
  {
    searchedAllElements = false;
    return true;
  }

  // look in every element, starting at elementNo-2
  element_no_t elementNoStart = (elementNo - 2 + nElements) % nElements;
  element_no_t elementNoEnd = (elementNo - 3 + nElements) % nElements;

  VLOG(3) << "elementNoStart: " << elementNoStart << ", elementNoEnd: " << elementNoEnd << ", nElements: " << nElements;
  searchedAllElements = true;

  // search among all elements
  for (element_no_t currentElementNo = elementNoStart; currentElementNo != elementNoEnd; currentElementNo++)
  {
    // restart with element 0
    if (currentElementNo == nElements)
    {
      currentElementNo = 0;
      if (elementNoEnd == currentElementNo)
        break;
    }

    if (this->pointIsInElement(point, currentElementNo, xi, residual, xiTolerance))
    {
      elementNo = currentElementNo;
      ghostMeshNo = -1;   // not a ghost mesh
      return true;
    }
  }
  return false;
}


template<int D,typename BasisFunctionType>
std::shared_ptr<FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> FunctionSpaceFindPosition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType,Mesh::UnstructuredDeformableOfDimension<D>>::
ghostMesh(Mesh::face_or_edge_t faceOrEdge)
{
  return nullptr;
}

} // namespace

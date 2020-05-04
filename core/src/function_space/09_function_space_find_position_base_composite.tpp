#include "function_space/09_function_space_find_position_base.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "mesh/face_t.h"

namespace FunctionSpace
{

// composite mesh
template<int D,typename BasisFunctionType>
bool FunctionSpaceStructuredFindPositionBase<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
findPosition(Vec3 point, element_no_t &elementNoLocal, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, bool startSearchInCurrentElement, double &residual, bool &searchedAllElements, double xiTolerance)
{
  const element_no_t nElements = this->nElementsLocal();
  VLOG(2) << "findPosition, elementNoLocal: " << elementNoLocal << ", ghostMeshNo: " << ghostMeshNo
    << ", startSearchInCurrentElement: " << startSearchInCurrentElement << ", xiTolerance: " << xiTolerance << ", xi: " << xi;

  // set starting no to 0 if it was not given and is thus arbitrarily initialized
  if (elementNoLocal < 0 || elementNoLocal >= nElements)
    elementNoLocal = 0;
  
  // get the submesh of the given elementNoLocal
  int subMeshOfElementSubmeshNo = 0;
  element_no_t elementOnMeshNoLocal = 0;
  this->meshPartition_->getSubMeshNoAndElementNoLocal(elementNoLocal, subMeshOfElementSubmeshNo, elementOnMeshNoLocal);

  bool firstIteration = true;
  searchedAllElements = false;

  // search among all elements in submeshes other than subMeshNo (because it was already searched there)
  for (int subMeshNo = subMeshOfElementSubmeshNo; subMeshNo != subMeshOfElementSubmeshNo || firstIteration; subMeshNo++)
  {
    if (subMeshNo == this->subFunctionSpaces_.size())
    {
      subMeshNo = 0;
      if (subMeshNo == subMeshOfElementSubmeshNo)
        break;
    }

    bool nodeFound = this->subFunctionSpaces_[subMeshNo]->findPosition(point, elementOnMeshNoLocal, ghostMeshNo, xi,
                                                                       startSearchInCurrentElement, residual, searchedAllElements, xiTolerance);
    if (nodeFound)
    {
      subMeshNoWherePointWasFound_ = subMeshNo;

      // transform elementOnMeshNoLocal, which is the element local no in the subMesh based numbering to the composite numbering
      elementNoLocal = this->meshPartition_->getElementNoLocalFromSubmesh(subMeshNo, elementOnMeshNoLocal);

      return true;
    }

    // if point was not found in first submesh, do not start from the given element no in the next submesh
    startSearchInCurrentElement = false;

    firstIteration = false;
  }

  VLOG(1) << "Could not find any containing element";
  return false;
}

} // namespace

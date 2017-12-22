#include "mesh/deformable.h"

#include <array>
#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include "easylogging++.h"

namespace Mesh
{

template<unsigned long D>
Deformable<D>::Deformable(PyObject *specificSettings) : Regular<D>(specificSettings)
{
  // compute number of nodes
  node_idx_t nNodes = this->nNodes();
  
  // create petsc vector that contains the node positions
  PetscErrorCode ierr;
  ierr = VecCreate(PETSC_COMM_WORLD, &nodePositions_);  CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) nodePositions_, "nodePositions"); CHKERRV(ierr);
  
  // initialize size of vector
  const int vectorSize = nNodes * 3;   // node Positions always contains 3 entries for every node
  ierr = VecSetSizes(nodePositions_, PETSC_DECIDE, vectorSize); CHKERRV(ierr);
  
  // set sparsity type and other options
  ierr = VecSetFromOptions(nodePositions_);  CHKERRV(ierr);
  
  // fill initial position from settings
  if (PythonUtility::containsKey(specificSettings, "nodePositions"))
  {
    int nodeDimension = PythonUtility::getOptionInt(specificSettings, "nodeDimension", 3, PythonUtility::ValidityCriterion::Between1And3);
    
    int inputVectorSize = nNodes * nodeDimension;
    std::vector<double> nodePositionsVector;
    PythonUtility::getOptionVector(specificSettings, "nodePositions", inputVectorSize, nodePositionsVector);
    
    LOG(DEBUG) << "nodeDimension: " << nodeDimension << ", expect input vector to have " << inputVectorSize << " entries.";

    // transform vector from (x,y) or (x) entries to (x,y,z) 
    if (nodeDimension < 3)
    {
      nodePositionsVector.resize(vectorSize);   // resize vector and value-initialize to 0
      for(int i=nNodes-1; i>=0; i--)
      {
        
        if (nodeDimension == 2)
          nodePositionsVector[i*3+1] = nodePositionsVector[i*nodeDimension+1];
        else
          nodePositionsVector[i*3+1] = 0;
        nodePositionsVector[i*3+0] = nodePositionsVector[i*nodeDimension+0];
        nodePositionsVector[i*3+2] = 0;
      }
    }
    
    PetscUtility::setVector(nodePositionsVector, nodePositions_);
  }
  else
  {
    // if node positions are not given in settings but physicalExtend, fill from that
    std::array<double, D> physicalExtend, meshWidth;
    physicalExtend = PythonUtility::getOptionArray<double, D>(specificSettings, "physicalExtend", 1.0, PythonUtility::Positive);
    
    for (unsigned int dimNo = 0; dimNo < D; dimNo++)
    {
      meshWidth[dimNo] = physicalExtend[dimNo] / this->nElements(dimNo);
      LOG(DEBUG) << "meshWidth["<<dimNo<<"] = "<<meshWidth[dimNo];
    }
    
    std::array<int, 3> indices;
    std::array<double, 3> position{0.,0.,0.};
    
    for (node_idx_t nodeNo = 0; nodeNo < nNodes; nodeNo++)
    {
      switch(D)
      {
      case 3:
        position[2] = meshWidth[2] * (int(nodeNo / (this->nNodes(0)*this->nNodes(1))));
      case 2:
        position[1] = meshWidth[1] * (int(nodeNo / this->nNodes(0)) % this->nNodes(1));
      case 1:
        position[0] = meshWidth[0] * (nodeNo % this->nNodes(0));
        
        break;
      }
      
      // set the indices where to store the position values
      indices[0] = nodeNo*3 + 0;
      indices[1] = nodeNo*3 + 1;
      indices[2] = nodeNo*3 + 2;
      
      VecSetValues(nodePositions_, 3, indices.data(), position.data(), INSERT_VALUES);
    }
  }
  
  // finish parallel assembly
  ierr = VecAssemblyBegin(nodePositions_); CHKERRV(ierr);
  ierr = VecAssemblyEnd(nodePositions_); CHKERRV(ierr);
  
  std::vector<double> nodePositions;
  PetscUtility::getVectorEntries(nodePositions_, nodePositions);

#ifndef NDEBUG  
  LOG(DEBUG) << "created nodePositions (" << nodePositions.size() << " entries): ";
  for (unsigned int i=0; i<nodePositions.size()/3; i++)
  {
    std::stringstream s;
    s << "(";
    for (unsigned int dimNo = 0; dimNo < 3; dimNo++)
    {
      if (dimNo != 0)
       s << ", ";
      s << nodePositions[i*3+dimNo];
    }
    s << ")";
    LOG(DEBUG) << s.str();
  }
#endif
  
}

template<unsigned long D>
void Deformable<D>::getGeometry(std::vector<double> &nodePositions)
{
  PetscUtility::getVectorEntries(nodePositions_, nodePositions);
}

template<unsigned long D>
Vec3 Deformable<D>::getGeometry(node_idx_t dofNo)
{
  Vec3 result;
  std::array<int,3> indices{dofNo*3+0, dofNo*3+1, dofNo*3+2};
  VecGetValues(nodePositions_, 3, indices.data(), result.data());
  return result;
}

};
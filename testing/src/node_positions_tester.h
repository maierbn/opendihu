#pragma once
  
#include <vector>
#include <map>
#include "utility/petsc_utility.h"
#include "mesh/mesh_manager.h"
#include "mesh/deformable.h"

namespace Mesh
{

class NodePositionsTester
{
public:
  static void compareNodePositions(
    DihuContext &dihuContext,
    std::string meshKey,
    std::vector<double> &referencePositions
  )
  {
    ASSERT_TRUE(dihuContext.meshManager()->meshes_.find(meshKey) != dihuContext.meshManager()->meshes_.end()) 
     << "Mesh with key \"" << meshKey << "\" was not found.";
    std::shared_ptr<Mesh> mesh = dihuContext.meshManager()->meshes_[meshKey];
    Vec nodePositions;
    
    if(mesh->dimension() == 1)
    {
      std::shared_ptr<Deformable<1>> deformableMesh = std::static_pointer_cast<Deformable<1>>(mesh);
      nodePositions = deformableMesh->nodePositions_;
    }
    else if(mesh->dimension() == 2)
    {
      std::shared_ptr<Deformable<2>> deformableMesh = std::static_pointer_cast<Deformable<2>>(mesh);
      nodePositions = deformableMesh->nodePositions_;
    }
    else if(mesh->dimension() == 3)
    {
      std::shared_ptr<Deformable<3>> deformableMesh = std::static_pointer_cast<Deformable<3>>(mesh);
      nodePositions = deformableMesh->nodePositions_;
    }
    std::vector<double> nodePositionsVector;
    PetscUtility::getVectorEntries(nodePositions, nodePositionsVector);
    
    ASSERT_EQ(nodePositionsVector.size(), referencePositions.size()) << "Wrong size of node positions vector.";
    
    for(unsigned int i=0; i<nodePositionsVector.size(); i++)
    {
      EXPECT_LE(fabs(nodePositionsVector[i]-referencePositions[i]), 1e-15) << "Node position " << i << " differs.";
    }
  }
};

};
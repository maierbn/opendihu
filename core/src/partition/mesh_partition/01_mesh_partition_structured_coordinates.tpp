#include "partition/mesh_partition/01_mesh_partition.h"

#include <cstdlib>
#include "utility/vector_operators.h"
#include "function_space/00_function_space_base_dim.h"

namespace Partition
{

template<typename MeshType,typename BasisFunctionType>
std::array<global_no_t,MeshType::dim()> MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getCoordinatesGlobal(node_no_t nodeNoLocal) const
{
  if (MeshType::dim() == 1)
  {
    return std::array<global_no_t,MeshType::dim()>({beginNodeGlobalNatural(0) + nodeNoLocal});
  }
  else if (MeshType::dim() == 2)
  {
    std::array<global_no_t,MeshType::dim()> coordinates;
    if (nodeNoLocal < nNodesLocalWithoutGhosts())
    {
      coordinates[0] = beginNodeGlobalNatural(0) + nodeNoLocal % nNodesLocalWithoutGhosts(0);
      coordinates[1] = beginNodeGlobalNatural(1) + nodeNoLocal / nNodesLocalWithoutGhosts(0);
    }
    else
    {
      if (hasFullNumberOfNodes(0))
      {
        if (hasFullNumberOfNodes(1))
        {
          // domain has no ghost nodes, should be handled by other if branch
          assert(false);
        }
        else
        {
          // domain has top ghost nodes
          global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts();
          coordinates[0] = beginNodeGlobalNatural(0) + ghostNo;
          coordinates[1] = beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1) - 1;
        }
      }
      else
      {
        if (hasFullNumberOfNodes(1))
        {
          // domain has right ghost nodes
          global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts();
          coordinates[0] = beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0) - 1;
          coordinates[1] = beginNodeGlobalNatural(1) + ghostNo;
        }
        else
        {
          // domain has right and top ghost nodes
          if (nodeNoLocal < nNodesLocalWithoutGhosts() + nNodesLocalWithoutGhosts(1))
          {
            // node is on right row of ghost nodes
            global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts();
            coordinates[0] = beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0) - 1;
            coordinates[1] = beginNodeGlobalNatural(1) + ghostNo;
          }
          else
          {
            // node is on top row of ghost nodes
            global_no_t ghostNo = nodeNoLocal - (nNodesLocalWithoutGhosts() + nNodesLocalWithGhosts(1) - 1);
            coordinates[0] = beginNodeGlobalNatural(0) + ghostNo;
            coordinates[1] = beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1) - 1;
          }
        }
      }
    }

    return coordinates;
  }
  else if (MeshType::dim() == 3)
  {
    std::array<global_no_t,MeshType::dim()> coordinates;
    if (nodeNoLocal < nNodesLocalWithoutGhosts())
    {
      coordinates[0] = beginNodeGlobalNatural(0) + nodeNoLocal % nNodesLocalWithoutGhosts(0);
      coordinates[1] = beginNodeGlobalNatural(1) + (nodeNoLocal % (nNodesLocalWithoutGhosts(0)*nNodesLocalWithoutGhosts(1))) / nNodesLocalWithoutGhosts(0);
      coordinates[2] = beginNodeGlobalNatural(2) + nodeNoLocal / (nNodesLocalWithoutGhosts(0)*nNodesLocalWithoutGhosts(1));
    }
    else
    {
      if (hasFullNumberOfNodes(0))
      {
        if (hasFullNumberOfNodes(1))
        {
          if (hasFullNumberOfNodes(2))
          {
            // domain has no ghost nodes, should be handled by other if branch
            assert(false);
          }
          else
          {
            // domain has only top (z+) ghost nodes
            int ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts();
            coordinates[0] = beginNodeGlobalNatural(0) + ghostNo % nNodesLocalWithGhosts(0);
            coordinates[1] = beginNodeGlobalNatural(1) + ghostNo / nNodesLocalWithGhosts(0);
            coordinates[2] = beginNodeGlobalNatural(2) + nNodesLocalWithGhosts(2) - 1;
          }
        }
        else
        {
          if (hasFullNumberOfNodes(2))
          {
            // domain has y+ ghost nodes
            global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts();
            coordinates[0] = beginNodeGlobalNatural(0) + ghostNo % nNodesLocalWithGhosts(0);
            coordinates[1] = beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1) - 1;
            coordinates[2] = beginNodeGlobalNatural(2) + ghostNo / nNodesLocalWithGhosts(0);
          }
          else
          {
            // domain has back (y+) and top (z+) ghost nodes
            if (nodeNoLocal < nNodesLocalWithoutGhosts() + nNodesLocalWithoutGhosts(0)*nNodesLocalWithoutGhosts(1))
            {
              // ghost is on back plane
              global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts();
              coordinates[0] = beginNodeGlobalNatural(0) + ghostNo % nNodesLocalWithoutGhosts(0);
              coordinates[1] = beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1) - 1;
              coordinates[2] = beginNodeGlobalNatural(2) + ghostNo / nNodesLocalWithoutGhosts(0);
            }
            else
            {
              // ghost is on top plane
              global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts() - nNodesLocalWithoutGhosts(0)*nNodesLocalWithoutGhosts(1);
              coordinates[0] = beginNodeGlobalNatural(0) + ghostNo % nNodesLocalWithoutGhosts(0);
              coordinates[1] = beginNodeGlobalNatural(1) + ghostNo / nNodesLocalWithoutGhosts(0);
              coordinates[2] = beginNodeGlobalNatural(2) + nNodesLocalWithGhosts(2) - 1;
            }
          }
        }
      }
      else
      {
        if (hasFullNumberOfNodes(1))
        {
          if (hasFullNumberOfNodes(2))
          {
            // domain has only right (x+) ghost nodes
            global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts();
            coordinates[0] = beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0) - 1;
            coordinates[1] = beginNodeGlobalNatural(1) + ghostNo % nNodesLocalWithoutGhosts(1);
            coordinates[2] = beginNodeGlobalNatural(2) + ghostNo / nNodesLocalWithoutGhosts(1);
          }
          else
          {
            // domain has right (x+) and top (z+) ghost nodes
            if (nodeNoLocal < nNodesLocalWithoutGhosts() + nNodesLocalWithoutGhosts(1)*nNodesLocalWithoutGhosts(2))
            {
              // ghost node is on right plane
              global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts();
              coordinates[0] = beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0) - 1;
              coordinates[1] = beginNodeGlobalNatural(1) + ghostNo % nNodesLocalWithoutGhosts(1);
              coordinates[2] = beginNodeGlobalNatural(2) + ghostNo / nNodesLocalWithoutGhosts(1);
            }
            else
            {
              // ghost node is on top plane
              global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts() - nNodesLocalWithoutGhosts(1)*nNodesLocalWithoutGhosts(2);
              coordinates[0] = beginNodeGlobalNatural(0) + ghostNo % nNodesLocalWithGhosts(0);
              coordinates[1] = beginNodeGlobalNatural(1) + ghostNo / nNodesLocalWithGhosts(0);
              coordinates[2] = beginNodeGlobalNatural(2) + nNodesLocalWithGhosts(2) - 1;
            }
          }
        }
        else
        {
          if (hasFullNumberOfNodes(2))
          {
            // domain has right (x+) and back (y+) ghost nodes
            global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts();
            global_no_t nNodesPerL = (nNodesLocalWithoutGhosts(1) + nNodesLocalWithGhosts(0));
            global_no_t ghostOnZPlaneNo = ghostNo % nNodesPerL;
            if (ghostOnZPlaneNo < nNodesLocalWithoutGhosts(1))
            {
              // ghost node is on right plane
              coordinates[0] = beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0) - 1;
              coordinates[1] = beginNodeGlobalNatural(1) + ghostOnZPlaneNo;
              coordinates[2] = beginNodeGlobalNatural(2) + ghostNo / nNodesPerL;
            }
            else
            {
              // ghost node is on back (y+) plane
              global_no_t ghostNoOnYPlane = ghostOnZPlaneNo - nNodesLocalWithoutGhosts(1);
              coordinates[0] = beginNodeGlobalNatural(0) + ghostNoOnYPlane;
              coordinates[1] = beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1) - 1;
              coordinates[2] = beginNodeGlobalNatural(2) + ghostNo / nNodesPerL;
            }
          }
          else
          {
            // domain has right (x+), back (y+) and top (z+) ghost nodes
            global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts();
            global_no_t nNodesPerL = (nNodesLocalWithoutGhosts(1) + nNodesLocalWithGhosts(0));
            if (ghostNo >= nNodesPerL * nNodesLocalWithoutGhosts(2))
            {
              // ghost node is on top plane
              ghostNo = ghostNo - nNodesPerL * nNodesLocalWithoutGhosts(2);
              coordinates[0] = beginNodeGlobalNatural(0) + ghostNo % nNodesLocalWithGhosts(0);
              coordinates[1] = beginNodeGlobalNatural(1) + ghostNo / nNodesLocalWithGhosts(0);
              coordinates[2] = beginNodeGlobalNatural(2) + nNodesLocalWithGhosts(2) - 1;
            }
            else
            {
              global_no_t ghostOnLNo = ghostNo % nNodesPerL;
              if (ghostOnLNo < nNodesLocalWithoutGhosts(1))
              {
                // ghost node is on right plane
                coordinates[0] = beginNodeGlobalNatural(0) + nNodesLocalWithGhosts(0) - 1;
                coordinates[1] = beginNodeGlobalNatural(1) + ghostOnLNo;
                coordinates[2] = beginNodeGlobalNatural(2) + ghostNo / nNodesPerL;
              }
              else
              {
                // ghost node is on back (y+) plane
                global_no_t ghostNoOnYPlane = ghostOnLNo - nNodesLocalWithoutGhosts(1);
                coordinates[0] = beginNodeGlobalNatural(0) + ghostNoOnYPlane;
                coordinates[1] = beginNodeGlobalNatural(1) + nNodesLocalWithGhosts(1) - 1;
                coordinates[2] = beginNodeGlobalNatural(2) + ghostNo / nNodesPerL;
              }
            }
          }
        }
      }
    }
    return coordinates;
  }
  else
  {
    assert(false);
  }
  return std::array<global_no_t,MeshType::dim()>({0});  // this never happens, but the intel compiler does not recognize it (gcc does)
}


template<typename MeshType,typename BasisFunctionType>
std::array<int,MeshType::dim()> MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getCoordinatesLocal(node_no_t nodeNoLocal) const
{
  std::array<global_no_t,MeshType::dim()> coordinatesGlobal = getCoordinatesGlobal(nodeNoLocal);
  std::array<int,MeshType::dim()> coordinatesLocal;

  if (MeshType::dim() == 1)
  {
    coordinatesLocal[0] = coordinatesGlobal[0] - beginNodeGlobalNatural(0);
  }
  else if (MeshType::dim() == 2)
  {
    coordinatesLocal[0] = coordinatesGlobal[0] - beginNodeGlobalNatural(0);
    coordinatesLocal[1] = coordinatesGlobal[1] - beginNodeGlobalNatural(1);
  }
  else if (MeshType::dim() == 3)
  {
    coordinatesLocal[0] = coordinatesGlobal[0] - beginNodeGlobalNatural(0);
    coordinatesLocal[1] = coordinatesGlobal[1] - beginNodeGlobalNatural(1);
    coordinatesLocal[2] = coordinatesGlobal[2] - beginNodeGlobalNatural(2);
  }
  else
  {
    assert(false);
  }
  return coordinatesLocal;
}

template<typename MeshType,typename BasisFunctionType>
std::array<int,MeshType::dim()> MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getCoordinatesLocal(std::array<global_no_t,MeshType::dim()> coordinatesGlobal, bool &isOnLocalDomain) const
{
  const int D = MeshType::dim();
  std::array<int,D> coordinatesLocal;
  isOnLocalDomain = true;
  for (int coordinateDirection = 0; coordinateDirection < D; coordinateDirection++)
  {
    coordinatesLocal[coordinateDirection] = coordinatesGlobal[coordinateDirection] - this->beginNodeGlobalNatural(coordinateDirection);

    // check if the computed local coordinate is in the local range of nodes
    if (coordinatesLocal[coordinateDirection] < 0 || coordinatesLocal[coordinateDirection] >= nNodesLocalWithoutGhosts(coordinateDirection))
    {
      isOnLocalDomain = false;
    }
  }

  return coordinatesLocal;

  //node_no_t nodeNoLocalFromCoordinates = functionSpace->getNodeNo(coordinatesLocal);
}

  //! get the local coordinates for a local element no
template<typename MeshType,typename BasisFunctionType>
std::array<int,MeshType::dim()> MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getElementCoordinatesLocal(element_no_t elementNoLocal) const
{
  std::array<int,MeshType::dim()> result;
  if (MeshType::dim() == 1)
  {
    result[0] = elementNoLocal;
  }
  else if (MeshType::dim() == 2)
  {
    result[1] = int(elementNoLocal / nElementsLocal_[0]);
    result[0] = elementNoLocal % nElementsLocal_[0];
  }
  else if (MeshType::dim() == 3)
  {
    result[2] = int(elementNoLocal / (nElementsLocal_[0] * nElementsLocal_[1]));
    result[1] = int((elementNoLocal % (nElementsLocal_[0] * nElementsLocal_[1])) / nElementsLocal_[0]);
    result[0] = elementNoLocal % nElementsLocal_[0];
  }
  else
  {
    assert(false);
  }
  return result;
}

}  // namespace

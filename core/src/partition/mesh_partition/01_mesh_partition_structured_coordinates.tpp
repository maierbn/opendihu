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
  if (nNodesLocalWithoutGhosts() == 0)
    return std::array<global_no_t,MeshType::dim()>({0});

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
            // for degenerate mesh partition
            coordinates[0] = 0;
            coordinates[1] = 0;
            coordinates[2] = 0;

            LOG(ERROR) << "degenerate mesh partition: " << nNodesLocalWithoutGhosts() << " local nodes";

            return coordinates;


            // domain has no ghost nodes, should be handled by other if branch
            //assert(false);
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
            if (nodeNoLocal < nNodesLocalWithoutGhosts() + nNodesLocalWithoutGhosts(0)*nNodesLocalWithoutGhosts(2))
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
              global_no_t ghostNo = nodeNoLocal - nNodesLocalWithoutGhosts() - nNodesLocalWithoutGhosts(0)*nNodesLocalWithoutGhosts(2);
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
  if (nNodesLocalWithoutGhosts() == 0)
    return std::array<int,MeshType::dim()>({0});

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

template<typename MeshType,typename BasisFunctionType>
std::array<global_no_t,MeshType::dim()> MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getCoordinatesGlobalOfNodeNoGlobalNatural(global_no_t nodeNoGlobalNatural) const
{
  std::array<global_no_t,MeshType::dim()> coordinatesGlobal;

  if (MeshType::dim() == 1)
  {
    coordinatesGlobal[0] = nodeNoGlobalNatural;
  }
  else if (MeshType::dim() == 2)
  {
    coordinatesGlobal[0] = nodeNoGlobalNatural % nNodesGlobal(0);
    coordinatesGlobal[1] = nodeNoGlobalNatural / nNodesGlobal(0);
  }
  else if (MeshType::dim() == 3)
  {
    coordinatesGlobal[0] = nodeNoGlobalNatural % nNodesGlobal(0);
    coordinatesGlobal[1] = (nodeNoGlobalNatural % (nNodesGlobal(0) * nNodesGlobal(1))) / nNodesGlobal(0);
    coordinatesGlobal[2] = nodeNoGlobalNatural / (nNodesGlobal(0) * nNodesGlobal(1));
  }

  return coordinatesGlobal;
}

//! get the local node no for its local coordinates, also works for ghost nodes
template<typename MeshType,typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>>::
getNodeNoLocal(std::array<int,MeshType::dim()> coordinatesLocal) const
{
  if (MeshType::dim() == 1)
  {
    return coordinatesLocal[0];
  }
  else if (MeshType::dim() == 2)
  {
    dof_no_t localX = coordinatesLocal[0];
    dof_no_t localY = coordinatesLocal[1];

    //VLOG(3) << "getNodeNo(" << coordinatesLocal << "), nNodesLocalWithoutGhosts: " << nNodesLocalWithoutGhosts()
    //  << ", nNodesLocalWithoutGhosts: (" << nNodesLocalWithoutGhosts(0) << "," << nNodesLocalWithoutGhosts(1) << ")";

    if (localX == nNodesLocalWithoutGhosts(0))   // point is on right ghost row
    {
      if (localY == nNodesLocalWithoutGhosts(1))  // point is on top ghost row
      {
        // top right ghost point
        //VLOG(3) << "  a: " << nNodesLocalWithGhosts()-1;
        return nNodesLocalWithGhosts()-1;
      }
      else
      {
        // on right ghost row
        //VLOG(3) << "  b: " << nNodesLocalWithoutGhosts() + localY;
        return nNodesLocalWithoutGhosts() + localY;
      }
    }
    else
    {
      if (localY == nNodesLocalWithoutGhosts(1))   // point is on top ghost row
      {
        // on top ghost row
        if (hasFullNumberOfNodes(0))
        {
          // there are only ghost on the top (y+)
          //VLOG(3) << "  c: " << nNodesLocalWithoutGhosts() + localX;
          return nNodesLocalWithoutGhosts() + localX;
        }
        else
        {
          // there are ghosts on the right (x+) and top (y+)
          //VLOG(3) << "  d: " << nNodesLocalWithoutGhosts() + nNodesLocalWithoutGhosts(1) + localX;
          return nNodesLocalWithoutGhosts() + nNodesLocalWithoutGhosts(1) + localX;
        }
      }
      else
      {
        // point is in interior
        //VLOG(3) << " e: " << nNodesLocalWithoutGhosts() + nNodesLocalWithoutGhosts(0)*localY + localX;
        return nNodesLocalWithoutGhosts(0)*localY + localX;
      }
    }
  }
  else if (MeshType::dim() == 3)
  {
    dof_no_t localX = coordinatesLocal[0];
    dof_no_t localY = coordinatesLocal[1];
    dof_no_t localZ = coordinatesLocal[2];

    /*VLOG(3) << "getNodeNo(" << coordinatesLocal << "), nNodesLocalWithoutGhosts: " << nNodesLocalWithoutGhosts()
      << ", nNodesLocalWithoutGhosts: (" << nNodesLocalWithoutGhosts(0) << "," << nNodesLocalWithoutGhosts(1) << "," << nNodesLocalWithoutGhosts(2) << ")"
      << ", hasFullNumberOfNodes: (" << hasFullNumberOfNodes(0) << "," << hasFullNumberOfNodes(1) << "," << hasFullNumberOfNodes(2) << ")";*/

    if (localZ == nNodesLocalWithoutGhosts(2))  // point is on top ghost row
    {
      if (hasFullNumberOfNodes(1)) // there is no ghosts at y+
      {
        if (hasFullNumberOfNodes(0)) // there is no ghosts at x+
        {
          // there are only ghosts at z+ and point is at z+
          return nNodesLocalWithoutGhosts()
            + nNodesLocalWithoutGhosts(0)*localY + localX;
        }
        else  // there are ghosts at x+
        {
          // there are ghosts at x+ and z+ and point is at z+
          return nNodesLocalWithoutGhosts()
            + nNodesLocalWithoutGhosts(1)*nNodesLocalWithoutGhosts(2)
            + nNodesLocalWithGhosts(0)*localY + localX;
        }
      }
      else  // there are ghosts at y+
      {
        if (hasFullNumberOfNodes(0)) // there is no ghosts at x+
        {
          // there are ghosts at y+ and z+ and point is at z+
          return nNodesLocalWithoutGhosts()
            + nNodesLocalWithoutGhosts(0)*nNodesLocalWithoutGhosts(2)
            + nNodesLocalWithoutGhosts(0)*localY + localX;
        }
        else  // there are ghosts at x+
        {
          // there are ghosts at x+, y+ and z+ and point is at z+
          const int nNodesPerL = nNodesLocalWithoutGhosts(0) + nNodesLocalWithGhosts(1);
          return nNodesLocalWithoutGhosts()
            + nNodesPerL*nNodesLocalWithoutGhosts(2)
            + nNodesLocalWithGhosts(0)*localY + localX;
        }
      }
    }
    else if (localX == nNodesLocalWithoutGhosts(0))   // point is on right ghost row (x+)
    {
      if (hasFullNumberOfNodes(1)) // there is no ghosts at y+
      {
        // there are only ghosts at x+
        return nNodesLocalWithoutGhosts()
            + nNodesLocalWithoutGhosts(1)*localZ + localY;
      }
      else  // there are ghosts at y+
      {
        const int nNodesPerL = nNodesLocalWithoutGhosts(0) + nNodesLocalWithGhosts(1);
        if (localY == nNodesLocalWithoutGhosts(1))  // point is on back ghost row (y+)
        {
          // ghost is at right back vertical edge
          return nNodesLocalWithoutGhosts()
            + nNodesPerL*localZ + nNodesLocalWithoutGhosts(1) + nNodesLocalWithoutGhosts(0);
        }
        else
        {
          // there are ghosts at x+ and y+
          return nNodesLocalWithoutGhosts()
              + nNodesPerL*localZ + localY;
        }
      }
    }
    else if (localY == nNodesLocalWithoutGhosts(1))  // point is on back ghost row (y+)
    {
      if (hasFullNumberOfNodes(0)) // there is no ghosts at x+
      {
        // there are only ghosts at y+
        return nNodesLocalWithoutGhosts()
            + nNodesLocalWithoutGhosts(0)*localZ + localX;
      }
      else  // there are ghosts at x+
      {
        // there are ghosts at x+ and y+
        const int nNodesPerL = nNodesLocalWithoutGhosts(0) + nNodesLocalWithGhosts(1);
        return nNodesLocalWithoutGhosts()
            + nNodesPerL*localZ + nNodesLocalWithoutGhosts(1) + localX;
      }
    }
    else
    {
      // ghost is in interior
      return nNodesLocalWithoutGhosts(0)*nNodesLocalWithoutGhosts(1)*localZ
        + nNodesLocalWithoutGhosts(0)*localY + localX;
    }
  }
#ifndef __PGI
  return 0;  // should not happen, but cray compiler does not recognize it
#endif
}

}  // namespace

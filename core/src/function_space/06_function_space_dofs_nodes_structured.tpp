#include "function_space/06_function_space_dofs_nodes_structured.h"

#include <Python.h>  // has to be the first included header
#include <cmath>
#include <array>

#include "easylogging++.h"

namespace FunctionSpace
{
  
template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
nNodesLocalWithGhosts(int coordinateDirection) const
{
  assert(this->meshPartition_);
  assert(this->meshPartition_->nElementsLocal(coordinateDirection) == this->nElementsPerCoordinateDirectionLocal(coordinateDirection));
 
  return this->meshPartition_->nNodesLocalWithGhosts(coordinateDirection);
}

template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
nNodesLocalWithoutGhosts(int coordinateDirection) const
{
  assert(this->meshPartition_);
  assert(this->meshPartition_->nElementsLocal(coordinateDirection) == this->nElementsPerCoordinateDirectionLocal(coordinateDirection));
 
  return this->meshPartition_->nNodesLocalWithoutGhosts(coordinateDirection);
}

template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
nNodesLocalWithGhosts() const
{
  assert(this->meshPartition_);
  return this->meshPartition_->nNodesLocalWithGhosts();
}

template<typename MeshType,typename BasisFunctionType>
node_no_t FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
nNodesLocalWithoutGhosts() const
{
  assert(this->meshPartition_);
  return this->meshPartition_->nNodesLocalWithoutGhosts();
}

template<typename MeshType,typename BasisFunctionType>
dof_no_t FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
nDofsLocalWithGhosts() const
{
  return nNodesLocalWithGhosts() * this->nDofsPerNode();
}

template<typename MeshType,typename BasisFunctionType>
dof_no_t FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
nDofsLocalWithoutGhosts() const
{
  return nNodesLocalWithoutGhosts() * this->nDofsPerNode();
}

template<typename MeshType,typename BasisFunctionType>
global_no_t FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
nNodesGlobal() const
{
  assert(this->meshPartition_);
  return this->meshPartition_->nNodesGlobal();
}

template<typename MeshType,typename BasisFunctionType>
global_no_t FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
nNodesGlobal(int coordinateDirection) const
{
  assert(this->meshPartition_);
  return this->meshPartition_->nNodesGlobal(coordinateDirection);
}

template<typename MeshType,typename BasisFunctionType>
global_no_t FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
nDofsGlobal() const
{
  return nNodesGlobal() * this->nDofsPerNode();
}

//! fill a vector with the node position entries, nodes will contain consecutively the (x,y,z) values of just all nodes, i.e. for Hermite not the derivatives
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
getNodePositions(std::vector<double> &nodes) const
{
  nodes.resize(this->nDofsLocalWithGhosts()*3);

  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < this->nDofsLocalWithGhosts(); nodeGlobalNo++)
  {

    node_no_t firstNodeDofGlobalNo = nodeGlobalNo*this->nDofsPerNode();

    int index = nodeGlobalNo*3;
    Vec3 position = this->geometryField_->getValue(firstNodeDofGlobalNo);
    nodes[index+0] = position[0];
    nodes[index+1] = position[1];
    nodes[index+2] = position[2];
  }
}

template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceDofsNodesStructured<MeshType,BasisFunctionType>::
setHermiteDerivatives()
{
  const int D = MeshType::dim();

  // for Hermite set derivatives as distances between nodes
  if (std::is_same<BasisFunctionType,BasisFunction::Hermite>::value)
  {
    // ensure that ghost buffers is filled with the correct values
    this->geometryField_->startGhostManipulation();

    // loop over nodes
    dof_no_t localDofNo = 0;
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < this->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      assert(localDofNo % this->nDofsPerNode() == 0);
      Vec3 nodePositionCurrent = this->geometryField_->getValue(localDofNo);

      // determine distance to neighbouring node in all directions (x-,x+,y-,y+,z-,z+)
      // distance to neighbours in x-direction
      node_no_t neighbour0MinusNodeNoLocal = this->getNeighbourNodeNoLocal(nodeNoLocal, Mesh::face_t::face0Minus);
      node_no_t neighbour0PlusNodeNoLocal = this->getNeighbourNodeNoLocal(nodeNoLocal, Mesh::face_t::face0Plus);

      Vec3 distance0Minus{0};
      Vec3 distance0Plus{0};
      Vec3 distance0{0};
      Vec3 distance1{0};
      Vec3 distance2{0};

      // determine distance to neighbouring node in 0- direction
      if (neighbour0MinusNodeNoLocal != -1)
      {
        // get dof No for node No, first dof of node (=0) for position value
        dof_no_t neighbour0MinusDofNoLocal = neighbour0MinusNodeNoLocal*this->nDofsPerNode() + 0;
        Vec3 nodePositionNeighbour0Minus = this->geometryField_->getValue(neighbour0MinusDofNoLocal);

        distance0Minus = nodePositionCurrent - nodePositionNeighbour0Minus;
      }

      // determine distance to neighbouring node in 0+ direction
      if (neighbour0PlusNodeNoLocal != -1)
      {
        // get dof No for node No, first dof of node (=0) for position value
        dof_no_t neighbour0PlusDofNoLocal = neighbour0PlusNodeNoLocal*this->nDofsPerNode() + 0;
        Vec3 nodePositionNeighbour0Plus = this->geometryField_->getValue(neighbour0PlusDofNoLocal);

        distance0Plus = nodePositionNeighbour0Plus - nodePositionCurrent;
      }

      VLOG(1) << "node local no. " << nodeNoLocal << ", \"0-\" " << neighbour0MinusNodeNoLocal << " (d: " << distance0Minus
        << "), \"0+\" " << neighbour0PlusNodeNoLocal << "(d: " << distance0Plus << ")";

      // average distance to both neighbours or take single value if there was only one neighbour found (because of end of domain)
      if (neighbour0MinusNodeNoLocal != -1 && neighbour0PlusNodeNoLocal != -1)
      {
        distance0 = 0.5*(distance0Minus + distance0Plus);
      }
      else if (neighbour0MinusNodeNoLocal != -1)
      {
        distance0 = distance0Minus;
      }
      else if (neighbour0PlusNodeNoLocal != -1)
      {
        distance0 = distance0Plus;
      }
      else
      {
        distance0 = Vec3{1.0};
      }

      VLOG(1) << "   final distance: " << distance0;

      if (D >= 2)
      {
        // distance to neighbours in y-direction
        node_no_t neighbour1MinusNodeNoLocal = this->getNeighbourNodeNoLocal(nodeNoLocal, Mesh::face_t::face1Minus);
        node_no_t neighbour1PlusNodeNoLocal = this->getNeighbourNodeNoLocal(nodeNoLocal, Mesh::face_t::face1Plus);

        Vec3 distance1Minus{0};
        Vec3 distance1Plus{0};

        // determine distance to neighbouring node in 1- direction
        if (neighbour1MinusNodeNoLocal != -1)
        {
          // get dof No for node No, first dof of node (=0) for position value
          dof_no_t neighbour1MinusDofNoLocal = neighbour1MinusNodeNoLocal*this->nDofsPerNode() + 0;
          Vec3 nodePositionNeighbour1Minus = this->geometryField_->getValue(neighbour1MinusDofNoLocal);

          distance1Minus = nodePositionCurrent - nodePositionNeighbour1Minus;
        }

        // determine distance to neighbouring node in 1+ direction
        if (neighbour1PlusNodeNoLocal != -1)
        {
          // get dof No for node No, first dof of node (=0) for position value
          dof_no_t neighbour1PlusDofNoLocal = neighbour1PlusNodeNoLocal*this->nDofsPerNode() + 0;
          Vec3 nodePositionNeighbour1Plus = this->geometryField_->getValue(neighbour1PlusDofNoLocal);

          distance1Plus = nodePositionNeighbour1Plus - nodePositionCurrent;
        }

        VLOG(1) << "node local no. " << nodeNoLocal << ", \"1-\" " << neighbour1MinusNodeNoLocal << " (d: " << distance1Minus
          << "), \"1+\" " << neighbour1PlusNodeNoLocal << "(d: " << distance1Plus << ")";

        // average distance to both neighbours or take single value if there was only one neighbour found (because of end of domain)
        if (neighbour1MinusNodeNoLocal != -1 && neighbour1PlusNodeNoLocal != -1)
        {
          distance1 = 0.5*(distance1Minus + distance1Plus);
        }
        else if (neighbour1MinusNodeNoLocal != -1)
        {
          distance1 = distance1Minus;
        }
        else if (neighbour1PlusNodeNoLocal != -1)
        {
          distance1 = distance1Plus;
        }
        else
        {
          distance1 = Vec3{1.0};
        }

        VLOG(1) << "   final distance: " << distance1;
      }

      if (D == 3)
      {
        // distance to neighbours in z-direction
        node_no_t neighbour2MinusNodeNoLocal = this->getNeighbourNodeNoLocal(nodeNoLocal, Mesh::face_t::face2Minus);
        node_no_t neighbour2PlusNodeNoLocal = this->getNeighbourNodeNoLocal(nodeNoLocal, Mesh::face_t::face2Plus);

        Vec3 distance2Minus{0};
        Vec3 distance2Plus{0};

        // determine distance to neighbouring node in 2- direction
        if (neighbour2MinusNodeNoLocal != -1)
        {
          // get dof No for node No, first dof of node (=0) for position value
          dof_no_t neighbour2MinusDofNoLocal = neighbour2MinusNodeNoLocal*this->nDofsPerNode() + 0;
          Vec3 nodePositionNeighbour2Minus = this->geometryField_->getValue(neighbour2MinusDofNoLocal);

          distance2Minus = nodePositionCurrent - nodePositionNeighbour2Minus;
        }

        // determine distance to neighbouring node in 2+ direction
        if (neighbour2PlusNodeNoLocal != -1)
        {
          // get dof No for node No, first dof of node (=0) for position value
          dof_no_t neighbour2PlusDofNoLocal = neighbour2PlusNodeNoLocal*this->nDofsPerNode() + 0;
          Vec3 nodePositionNeighbour2Plus = this->geometryField_->getValue(neighbour2PlusDofNoLocal);

          distance2Plus = nodePositionNeighbour2Plus - nodePositionCurrent;
        }

        VLOG(1) << "node local no. " << nodeNoLocal << ", \"2-\" " << neighbour2MinusNodeNoLocal << " (d: " << distance2Minus
          << "), \"2+\" " << neighbour2PlusNodeNoLocal << "(d: " << distance2Plus << ")";

        // average distance to both neighbours or take single value if there was only one neighbour found (because of end of domain)
        if (neighbour2MinusNodeNoLocal != -1 && neighbour2PlusNodeNoLocal != -1)
        {
          distance2 = 0.5*(distance2Minus + distance2Plus);
        }
        else if (neighbour2MinusNodeNoLocal != -1)
        {
          distance2 = distance2Minus;
        }
        else if (neighbour2PlusNodeNoLocal != -1)
        {
          distance2 = distance2Plus;
        }
        else
        {
          distance2 = Vec3{1.0};
        }

        VLOG(1) << "   final distance: " << distance2;
      }

      // advance localDof No to first non-positional (derivative) value
      localDofNo++;

      if (D == 1)
      {
        this->geometryField_->setValue(localDofNo, distance0);
        localDofNo++;
      }
      else if (D == 2)
      {
        // y x
        this->geometryField_->setValue(localDofNo, 1*distance0);
        localDofNo++;

        this->geometryField_->setValue(localDofNo, distance1*1);
        localDofNo++;

        this->geometryField_->setValue(localDofNo, distance1*distance0);
        localDofNo++;
      }
      else if (D == 3)
      {
        // z y x
        this->geometryField_->setValue(localDofNo, 1*1*distance0);
        localDofNo++;

        this->geometryField_->setValue(localDofNo, 1*distance1*1);
        localDofNo++;

        this->geometryField_->setValue(localDofNo, 1*distance1*distance0);
        localDofNo++;

        this->geometryField_->setValue(localDofNo, distance2*1*1);
        localDofNo++;

        this->geometryField_->setValue(localDofNo, distance2*1*distance0);
        localDofNo++;

        this->geometryField_->setValue(localDofNo, distance2*distance1*1);
        localDofNo++;

        this->geometryField_->setValue(localDofNo, distance2*distance1*distance0);
        localDofNo++;
      }
    }
  }
}

} // namespace

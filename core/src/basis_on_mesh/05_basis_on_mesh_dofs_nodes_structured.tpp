#include "basis_on_mesh/05_basis_on_mesh_dofs_nodes.h"

#include <Python.h>  // has to be the first included header
#include <cmath>
#include <array>

#include "easylogging++.h"
#include "field_variable/field_variable.h"
#include "utility/petsc_utility.h"
#include "basis_on_mesh/00_basis_on_mesh_base_dim.h"

namespace BasisOnMesh
{

// constructor
template<int D,typename BasisFunctionType>
BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
BasisOnMeshDofsNodes(PyObject *specificSettings, bool noGeometryField) :
  BasisOnMeshGeometry<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>(specificSettings)
{
  LOG(DEBUG) << "constructor BasisOnMeshDofsNodes StructuredDeformable, noGeometryField_="<<this->noGeometryField_;

  this->noGeometryField_ = noGeometryField;
  std::vector<double> nodePositions;
  if (!this->noGeometryField_)
  {
    this->parseNodePositionsFromSettings(specificSettings, nodePositions);
    this->setGeometryField(nodePositions);
  }
}

template<int D,typename BasisFunctionType>
BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
BasisOnMeshDofsNodes(const std::vector<Vec3> &nodePositions) :
  BasisOnMeshGeometry<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>(NULL)
{
  this->noGeometryField_ = false;

  std::vector<double> sequentialNodePositions;   // node positions in a scalar vector, as needed by setGeometryField
  sequentialNodePositions.reserve(nodePositions.size() * D);

  for (Vec3 &vector : nodePositions)
  {
    for (int i = 0; i < 3; i++)
      sequentialNodePositions.push_back(vector[i]);
  }
  this->setGeometryField(sequentialNodePositions);
}


// read in config nodes
template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
parseNodePositionsFromSettings(PyObject *specificSettings, std::vector<double> &nodePositions)
{
  // compute number of nodes
  node_no_t nNodes = this->nNodes();

  const int vectorSize = nNodes*3;
  nodePositions.resize(vectorSize);   // resize vector and value-initialize to 0

  // fill initial position from settings
  if (PythonUtility::hasKey(specificSettings, "nodePositions"))
  {
    bool nodesStoredAsLists = false;

    // check if the node positions are stored as list, e.g. [[x,y,z], [x,y,z],...]
    PyObject *nodePositionsListPy = PythonUtility::getOptionPyObject(specificSettings, "nodePositions");
    if (PyList_Check(nodePositionsListPy))
    {
      if (PyList_Size(nodePositionsListPy) > 0)
      {
        PyObject *firstItemPy = PyList_GetItem(nodePositionsListPy, (Py_ssize_t)0);
        if (PyList_Check(firstItemPy))
        {
          nodesStoredAsLists = true;
        }
      }
    }

    LOG(DEBUG) << "nodePositions: " << nodePositionsListPy << ", nodesStoredAsLists=" << nodesStoredAsLists;

    if (nodesStoredAsLists)
    {
      node_no_t nNodesInList = PyList_Size(nodePositionsListPy);
      node_no_t nodeNo = 0;
      for (; nodeNo < nNodesInList; nodeNo++)
      {
        // extract single node position, e.g. [x,y]
        PyObject *itemNodePositionPy = PyList_GetItem(nodePositionsListPy, (Py_ssize_t)nodeNo);
        if (PyList_Check(itemNodePositionPy))
        {
          // parse components of node, e.g. x
          int i = 0;
          for (; i < std::min(3,(int)PyList_Size(itemNodePositionPy)); i++)
          {
            PyObject *pointComponentPy = PyList_GetItem(itemNodePositionPy, (Py_ssize_t)i);
            nodePositions[3*nodeNo + i] = PythonUtility::convertFromPython<double>(pointComponentPy, 0.0);
          }

          // set the rest of the values that were not specified to 0.0, e.g. z=0.0
          for (; i < 3; i++)
          {
            nodePositions[3*nodeNo + i] = 0.0;
          }

          LOG(DEBUG) << "(1) set node " << nodeNo << "[" << nodePositions[3*nodeNo + 0] << "," << nodePositions[3*nodeNo + 1] << "," << nodePositions[3*nodeNo + 2] << "]";
        }
        else
        {
          // if the entry is not a list like [x,y,z] but a single value, assume it is the x value
          double value = PythonUtility::convertFromPython<double>(itemNodePositionPy, 0.0);
          nodePositions[3*nodeNo + 0] = value;
          nodePositions[3*nodeNo + 1] = 0.0;
          nodePositions[3*nodeNo + 2] = 0.0;

          LOG(DEBUG) << "(2) set node " << nodeNo << "[" << nodePositions[3*nodeNo + 0] << "," << nodePositions[3*nodeNo + 1] << "," << nodePositions[3*nodeNo + 2] << "]";
        }
      }

      if (nodeNo < nNodes)
      {
        LOG(WARNING) << "Expected " << nNodes << " nodes, nodePositions contains only " << nodeNo << " nodes (" << nNodes - nodeNo << " missing).";
      }

      // fill rest of values with 0,0,0
      for (; nodeNo < nNodes; nodeNo++)
      {
        nodePositions[3*nodeNo + 0] = 0.0;
        nodePositions[3*nodeNo + 1] = 0.0;
        nodePositions[3*nodeNo + 2] = 0.0;

        LOG(DEBUG) << "(3) set node " << nodeNo << "[" << nodePositions[3*nodeNo + 0] << "," << nodePositions[3*nodeNo + 1] << "," << nodePositions[3*nodeNo + 2] << "]";
      }
    }
    else
    {
      // nodes are stored as contiguous array, e.g. [x,y,z,x,y,z] or [x,y,x,y,x,y,...]

      int nodeDimension = PythonUtility::getOptionInt(specificSettings, "nodeDimension", 3, PythonUtility::ValidityCriterion::Between1And3);

      int inputVectorSize = nNodes * nodeDimension;
      PythonUtility::getOptionVector(specificSettings, "nodePositions", inputVectorSize, nodePositions);

      LOG(DEBUG) << "nodeDimension: " << nodeDimension << ", expect input vector to have " << nNodes << "*" << nodeDimension << "=" << inputVectorSize << " entries.";

      // transform vector from (x,y) or (x) entries to (x,y,z)
      if (nodeDimension < 3)
      {
        nodePositions.resize(vectorSize);   // resize vector and value-initialize to 0
        for(int i=nNodes-1; i>=0; i--)
        {

          if (nodeDimension == 2)
            nodePositions[i*3+1] = nodePositions[i*nodeDimension+1];
          else
            nodePositions[i*3+1] = 0;
          nodePositions[i*3+0] = nodePositions[i*nodeDimension+0];
          nodePositions[i*3+2] = 0;
        }
      }
    }
  }
  else   // there was no "nodePositions" given in config, use physicalExtent instead
  {
    // if node positions are not given in settings but physicalExtent, fill from that
    std::array<double, D> physicalExtent, meshWidth;
    physicalExtent = PythonUtility::getOptionArray<double, D>(specificSettings, "physicalExtent", 1.0, PythonUtility::Positive);

    for (unsigned int dimNo = 0; dimNo < D; dimNo++)
    {
      meshWidth[dimNo] = physicalExtent[dimNo] /
        (this->nElementsPerCoordinateDirection(dimNo) * (BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement()-1));
      LOG(DEBUG) << "meshWidth["<<dimNo<<"] = "<<meshWidth[dimNo];
    }

    std::array<double, 3> position{0.,0.,0.};

    LOG(DEBUG) << "nNodes: " << nNodes << ", vectorSize: " << vectorSize;

    for (node_no_t nodeNo = 0; nodeNo < nNodes; nodeNo++)
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


      // store the position values in nodePositions
      for (int i=0; i<3; i++)
        nodePositions[nodeNo*3 + i] = position[i];
    }
  }

}

// create geometry field from config nodes
template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
setGeometryField(std::vector<double> &nodePositions)
{
  LOG(DEBUG) << " BasisOnMesh StructuredDeformable, setGeometryField";

  // compute number of dofs
  dof_no_t nDofs = this->nDofs();

  // create petsc vector that contains the node positions
  Vec values;
  PetscErrorCode ierr;
  ierr = VecCreate(PETSC_COMM_WORLD, &values);  CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) values, "geometry"); CHKERRV(ierr);

  // initialize size of vector
  const int vectorSize = nDofs * 3;   // dofs always contain 3 entries for every entry (x,y,z)
  ierr = VecSetSizes(values, PETSC_DECIDE, vectorSize); CHKERRV(ierr);

  // set sparsity type and other options
  ierr = VecSetFromOptions(values);  CHKERRV(ierr);

  // fill geometry vector from nodePositions, initialize non-node position entries to 0 (for Hermite)
  std::vector<double> geometryValues(vectorSize, 0.0);

  int geometryValuesIndex = 0;
  int nodePositionsIndex = 0;
  // loop over nodes
  for (node_no_t nodeNo = 0; nodeNo < this->nNodes(); nodeNo++)
  {
    // assign node position as first dof of the node
    geometryValues[geometryValuesIndex+0] = nodePositions[nodePositionsIndex+0];
    geometryValues[geometryValuesIndex+1] = nodePositions[nodePositionsIndex+1];
    geometryValues[geometryValuesIndex+2] = nodePositions[nodePositionsIndex+2];
    geometryValuesIndex += 3;
    nodePositionsIndex += 3;

    // set entries to 0 for rest of dofs at this node
    for (int dofIndex = 1; dofIndex < this->nDofsPerNode(); dofIndex++)
    {
      geometryValues[geometryValuesIndex+0] = 0;
      geometryValues[geometryValuesIndex+1] = 0;
      geometryValues[geometryValuesIndex+2] = 0;
      geometryValuesIndex += 3;
    }
  }

  LOG(DEBUG) << "setGeometryField, geometryValues: " << geometryValues;

  PetscUtility::setVector(geometryValues, values);

  // finish parallel assembly
  ierr = VecAssemblyBegin(values); CHKERRV(ierr);
  ierr = VecAssemblyEnd(values); CHKERRV(ierr);

  bool isGeometryField = true;   // if the field is a geometry field
  // set geometry field
  this->geometryField_ = std::make_unique<GeometryFieldType>();
  std::vector<std::string> componentNames{"x", "y", "z"};
  int nEntries = nDofs * 3;   // 3 components (x,y,z) per dof
  this->geometryField_->set("geometry", componentNames, this->nElementsPerCoordinateDirection_, nEntries, isGeometryField, values);
}

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodes() const
{
  int result = 1;
  for (int i=0; i<D; i++)
    result *= nNodes(i);
  return result;
}

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nNodes(int dimension) const
{
  return this->nElementsPerCoordinateDirection(dimension) * BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement() + 1;
}

template<int D,typename BasisFunctionType>
dof_no_t BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
nDofs() const
{
  return nNodes() * this->nDofsPerNode();
}


template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
getNodePositions(std::vector<double> &nodes) const
{
  nodes.resize(this->nNodes()*3);

  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < this->nNodes(); nodeGlobalNo++)
  {

    node_no_t firstNodeDofGlobalNo = nodeGlobalNo*this->nDofsPerNode();

    int index = nodeGlobalNo*3;
    Vec3 position = this->geometryField_->getValue(firstNodeDofGlobalNo);
    nodes[index+0] = position[0];
    nodes[index+1] = position[1];
    nodes[index+2] = position[2];
  }
}


};  // namespace
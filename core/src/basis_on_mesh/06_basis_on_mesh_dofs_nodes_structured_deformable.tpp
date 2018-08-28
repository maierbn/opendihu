#include "basis_on_mesh/06_basis_on_mesh_dofs_nodes.h"

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
BasisOnMeshDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, PyObject *specificSettings, bool noGeometryField) :
  BasisOnMeshDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>(partitionManager, specificSettings)
{
  LOG(DEBUG) << "constructor BasisOnMeshDofsNodes StructuredDeformable, noGeometryField_=" <<this->noGeometryField_;

  this->noGeometryField_ = noGeometryField;
}

template<int D,typename BasisFunctionType>
BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
BasisOnMeshDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, const std::vector<Vec3> &localNodePositions, const std::array<element_no_t,D> nElementsPerCoordinateDirection) :
  BasisOnMeshDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>(NULL)
{
  LOG(DEBUG) << "constructor BasisOnMeshDofsNodes StructuredDeformable, from " << localNodePositions.size() << " localNodePositions";
 
  this->noGeometryField_ = false;
  this->nElementsPerCoordinateDirectionLocal_ = nElementsPerCoordinateDirection;
  LOG(DEBUG) << "set number of elements per coordinate direction: " << this->nElementsPerCoordinateDirectionLocal_;

  localNodePositions_.reserve(localNodePositions.size() * D);

  for (const Vec3 &vector : localNodePositions)
  {
    for (int i = 0; i < 3; i++)
      localNodePositions_.push_back(vector[i]);
  }
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
initialize()
{ 
  // create meshPartition and redistribute elements if necessary, this needs information about mesh size
  BasisOnMeshPartition<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::initialize();
 
  // if this mesh does not have a geometry field, do nothing further
  if (this->noGeometryField_)
    return;

  // setup geometry field
  // parse node positions from python config
  this->parseNodePositionsFromSettings(this->specificSettings_);
 
  // pass a shared "this" pointer to the geometryField 
  // retrieve "this" pointer and convert to downwards pointer of most derived class "BasisOnMesh"
  std::shared_ptr<BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> thisMesh
    = std::static_pointer_cast<BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>(this->shared_from_this());

  assert(thisMesh != nullptr);
  
  // create empty field variable for geometry field
  std::vector<std::string> componentNames{"x", "y", "z"};
  this->geometryField_ = std::make_unique<GeometryFieldType>(thisMesh, "geometry", componentNames);
  
  // assign values of geometry field
  this->setGeometryFieldValues();
}

// read in config nodes
template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
parseNodePositionsFromSettings(PyObject *specificSettings)
{
  LOG(TRACE) << "BasisOnMeshDofsNodes<structuredDeformable> parseNodePositionsFromSettings";
  
  // compute number of nodes
  node_no_t nNodesLocal = this->nNodesLocalWithoutGhosts();
  global_no_t nNodesGlobal = this->nNodesGlobal();
  global_no_t nNodes;

  // if the given information about the mesh in config is for the global mesh
  bool inputMeshIsGlobal = PythonUtility::getOptionBool(specificSettings, "inputMeshIsGlobal", true);
  LOG(DEBUG) << "inputMeshIsGlobal: " << std::boolalpha << inputMeshIsGlobal;
  
  global_no_t vectorSize;
  
  if (inputMeshIsGlobal)
  {
    vectorSize = nNodesGlobal*3;
    nNodes = nNodesGlobal;
  }
  else 
  {
    vectorSize = nNodesLocal*3;
    nNodes = nNodesLocal;
  }
  
  VLOG(1) << "nNodesLocal: " << nNodesLocal << ", nNodesGlobal: " << nNodesGlobal << ", vectorSize: " << vectorSize;
  
  localNodePositions_.resize(vectorSize);   // resize vector and value-initialize to 0
  // The vector localNodePositions_ is filled with all available node positions in this method, for inputMeshIsGlobal these are the global values of the global domain.
  // At the end of this method, non-local entries are removed if inputMeshIsGlobal.

  
  // fill initial position from settings
  if (PythonUtility::hasKey(specificSettings, "nodePositions"))
  {
    bool nodesStoredAsLists = false;
    VLOG(1) << "specificSettings has \"nodePositions\"";

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

    VLOG(1) << "nodePositions: " << nodePositionsListPy << ", nodesStoredAsLists=" << nodesStoredAsLists;

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
            localNodePositions_[3*nodeNo + i] = PythonUtility::convertFromPython<double>(pointComponentPy, 0.0);
          }

          // set the rest of the values that were not specified to 0.0, e.g. z=0.0
          for (; i < 3; i++)
          {
            localNodePositions_[3*nodeNo + i] = 0.0;
          }

          VLOG(2) << "(1) set node " << nodeNo << "[" << localNodePositions_[3*nodeNo + 0] << "," << localNodePositions_[3*nodeNo + 1] << "," << localNodePositions_[3*nodeNo + 2] << "]";
        }
        else
        {
          // if the entry is not a list like [x,y,z] but a single value, assume it is the x value
          double value = PythonUtility::convertFromPython<double>(itemNodePositionPy, 0.0);
          localNodePositions_[3*nodeNo + 0] = value;
          localNodePositions_[3*nodeNo + 1] = 0.0;
          localNodePositions_[3*nodeNo + 2] = 0.0;

          VLOG(2) << "(2) set node " << nodeNo << "[" << localNodePositions_[3*nodeNo + 0] << "," << localNodePositions_[3*nodeNo + 1] << "," << localNodePositions_[3*nodeNo + 2] << "]";
        }
      }

      if (nodeNo < nNodes)
      {
        LOG(WARNING) << "Expected " << nNodes << " nodes, localNodePositions_ contains only " << nodeNo << " nodes (" << nNodes - nodeNo << " missing).";
      }

      // fill rest of values with 0,0,0
      for (; nodeNo < nNodes; nodeNo++)
      {
        localNodePositions_[3*nodeNo + 0] = 0.0;
        localNodePositions_[3*nodeNo + 1] = 0.0;
        localNodePositions_[3*nodeNo + 2] = 0.0;

        VLOG(2) << "(3) set node " << nodeNo << "[" << localNodePositions_[3*nodeNo + 0] << "," << localNodePositions_[3*nodeNo + 1] << "," << localNodePositions_[3*nodeNo + 2] << "]";
      }
    }
    else
    {
      // nodes are stored as contiguous array, e.g. [x,y,z,x,y,z] or [x,y,x,y,x,y,...]

      int nodeDimension = PythonUtility::getOptionInt(specificSettings, "nodeDimension", 3, PythonUtility::ValidityCriterion::Between1And3);

      int inputVectorSize = nNodesLocal * nodeDimension;
      PythonUtility::getOptionVector(specificSettings, "nodePositions", inputVectorSize, localNodePositions_);

      LOG(DEBUG) << "nodeDimension: " << nodeDimension << ", expect input vector to have " << nNodesLocal << "*" << nodeDimension << "=" << inputVectorSize << " entries.";

      // transform vector from (x,y) or (x) entries to (x,y,z)
      if (nodeDimension < 3)
      {
        localNodePositions_.resize(vectorSize);   // resize vector and value-initialize to 0
        for(int i = nNodesLocal-1; i >= 0; i--)
        {

          if (nodeDimension == 2)
          {
            localNodePositions_[i*3+1] = localNodePositions_[i*nodeDimension+1];
          }
          else
          {
            localNodePositions_[i*3+1] = 0;
          }
          localNodePositions_[i*3+0] = localNodePositions_[i*nodeDimension+0];
          localNodePositions_[i*3+2] = 0;
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
      // depending on inputMeshIsGlobal interpret physicalExtent as global or local physical extent
      if (inputMeshIsGlobal)
      {
        meshWidth[dimNo] = physicalExtent[dimNo] /
          (this->nElementsPerCoordinateDirectionGlobal(dimNo) * (BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement()-1));
      }
      else 
      {
        meshWidth[dimNo] = physicalExtent[dimNo] /
          (this->nElementsPerCoordinateDirectionLocal(dimNo) * (BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement()-1));
      }
      LOG(DEBUG) << "meshWidth[" <<dimNo<< "] = " <<meshWidth[dimNo];
    }
   
    VLOG(1) << "specificSettings has no \"nodePositions\", use physicalExtent: " << physicalExtent << ", meshWidth: " << meshWidth;
    
    std::array<double, 3> position{0.,0.,0.};

    // compute absolute node positions
    global_no_t nNodesInXDirection;
    global_no_t nNodesInYDirection;
    global_no_t nNodesInZDirection;
    
    double offsetX;
    double offsetY;
    double offsetZ;
    
    // initialize variables
    if (inputMeshIsGlobal)
    {
      switch(D)
      {
      case 3:
        nNodesInZDirection = this->nNodesGlobal(2);
      case 2:
        nNodesInYDirection = this->nNodesGlobal(1);
      case 1:
        nNodesInXDirection = this->nNodesGlobal(0);
      }
      
      offsetX = 0;
      offsetY = 0;
      offsetZ = 0;
    }
    else 
    {
      switch(D)
      {
      case 3:
        nNodesInZDirection = this->nNodesLocalWithoutGhosts(2);
        offsetZ = this->meshPartition_->beginNodeGlobal(2);
      case 2:
        nNodesInYDirection = this->nNodesLocalWithoutGhosts(1);
        offsetY = this->meshPartition_->beginNodeGlobal(1);
      case 1:
        nNodesInXDirection = this->nNodesLocalWithoutGhosts(0);
        offsetX = this->meshPartition_->beginNodeGlobal(0);
      }
    }
    
    global_no_t nodeX, nodeY, nodeZ;
    LOG(DEBUG) << "nNodes: " << nNodes << ", global: " << nNodesGlobal << ", local: " << nNodesLocal << ", vectorSize: " << vectorSize;

    for (global_no_t nodeNo = 0; nodeNo < nNodes; nodeNo++)
    {
      switch(D)
      {
      case 3:
        nodeZ = global_no_t(nodeNo / (nNodesInXDirection*nNodesInYDirection));
        position[2] = offsetZ + meshWidth[2] * nodeZ;
      case 2:
        nodeY = global_no_t(nodeNo / nNodesInXDirection) % nNodesInYDirection;
        position[1] = offsetY + meshWidth[1] * nodeY;
      case 1:
        nodeX = nodeNo % nNodesInXDirection;
        position[0] = offsetX + meshWidth[0] * nodeX;

        break;
      }

      // store the position values in nodePositions
      for (int i = 0; i < 3; i++)
      {
        localNodePositions_[nodeNo*3 + i] = position[i];
      }
    }
  }

  // if parsed node positions in vector localNodePositions_ actually contains global node positions, extract local positions
  if (inputMeshIsGlobal)
  {
    this->meshPartition_->extractLocalNodes(localNodePositions_);
  }
}

// create geometry field from config nodes
template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodes<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::
setGeometryFieldValues()
{
  LOG(DEBUG) << " BasisOnMesh StructuredDeformable, setGeometryField, size of nodePositions vector: " << localNodePositions_.size();

  // compute number of (local) dofs
  dof_no_t nDofsLocal = this->nDofsLocalWithoutGhosts();
  
  // fill geometry vector from nodePositions, initialize non-node position entries to 0 (for Hermite)
  std::vector<Vec3> geometryValues(nDofsLocal, Vec3{0.0});

  int geometryValuesIndex = 0;
  int nodePositionsIndex = 0;
  // loop over nodes
  for (node_no_t nodeNo = 0; nodeNo < this->nNodesLocalWithoutGhosts(); nodeNo++)
  {
    // assign node position as first dof of the node
    geometryValues[geometryValuesIndex] = Vec3
    { 
      localNodePositions_[nodePositionsIndex+0], 
      localNodePositions_[nodePositionsIndex+1], 
      localNodePositions_[nodePositionsIndex+2]
    };
    geometryValuesIndex++;
    nodePositionsIndex += 3;

    // set entries to 0 for rest of dofs at this node (versions)
    for (int dofIndex = 1; dofIndex < this->nDofsPerNode(); dofIndex++)
    {
      geometryValuesIndex++;
    }
  }

  LOG(DEBUG) << "setGeometryField, geometryValues: " << geometryValues.size();

  // set values for node positions as geometry field 
  this->geometryField_->setValuesWithoutGhosts(geometryValues);
  this->geometryField_->finishVectorManipulation();
}


};  // namespace

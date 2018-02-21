#include "basis_on_mesh/03_basis_on_mesh_dofs.h"

#include <cmath>
#include <array>
#include <string>
#include <map>
#include <cassert>

#include "easylogging++.h"
#include "utility/string_utility.h"
#include "utility/math_utility.h"
#include "field_variable/field_variable_regular_fixed.h"
#include "field_variable/field_variable_structured_deformable.h"
#include "field_variable/field_variable_unstructured_deformable.h"

namespace BasisOnMesh
{
 
using namespace StringUtility;

// element-local dofIndex to global dofNo for 1D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  return BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>>::getDofNo(this->nElementsPerCoordinateDirection_, elementNo, dofIndex);
}

template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getDofNo(std::array<element_no_t, MeshType::dim()> nElementsPerCoordinateDirection, element_no_t elementNo, int dofIndex)
{
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0,1 2,3
  // averageNDofsPerElement: 
  // 1         2            2
  return BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNDofsPerElement() * elementNo + dofIndex;
}

//! get all dofs of a specific node for 1D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<int> dofGlobalNos) const
{
  for (int i=0; i<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos.push_back(BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i);
  }
}

// element-local dofIndex to global dofNo for 2D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  return BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>>::getDofNo(this->nElementsPerCoordinateDirection_, elementNo, dofIndex);
}

// element-local dofIndex to global dofNo for 2D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getDofNo(std::array<element_no_t, MeshType::dim()> nElements, element_no_t elementNo, int dofIndex)
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      45 67
  // 0 1       0 1 2      01 23
  // nDofsPerBasis:
  // 2         3          4
  // averageNDofsPerElement:
  // 1         4          2
  
  int averageNDofsPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNDofsPerElement();
  int dofsPerRow = (averageNDofsPerElement1D * nElements[0] + BasisFunctionType::nDofsPerNode());
  int elementX = int(elementNo % nElements[0]);
  int elementY = int(elementNo / nElements[0]);
  int localX = dofIndex % BasisFunctionType::nDofsPerBasis();
  int localY = int(dofIndex / BasisFunctionType::nDofsPerBasis());
  
  VLOG(2) << "  dof " << elementNo << ":" << dofIndex << ", element: ("<<elementX<<","<<elementY<<"), dofsPerRow="<<dofsPerRow<<", local: ("<<localX<<","<<localY<<")";
  
  return dofsPerRow * (elementY * averageNDofsPerElement1D + localY) 
    + averageNDofsPerElement1D * elementX + localX;
}

//! get all dofs of a specific node for 2D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<int> dofGlobalNos) const
{
  for (int i=0; i<BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos.push_back(BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i);
  }
}
  
// element-local dofIndex to global dofNo for 3D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  return BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>>::getDofNo(this->nElementsPerCoordinateDirection_, elementNo, dofIndex);
}
 
// element-local dofIndex to global dofNo for 3D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getDofNo(std::array<element_no_t, MeshType::dim()> nElements, element_no_t elementNo, int dofIndex)
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      45 67
  // 0 1       0 1 2      01 23
  // nDofsPerBasis:
  // 2         3          4
  // averageNDofsPerElement:
  // 1         4          2
 
  int averageNDofsPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNDofsPerElement();
  int dofsPerRow = (averageNDofsPerElement1D * nElements[0] + BasisFunctionType::nDofsPerNode());
  int dofsPerPlane = (averageNDofsPerElement1D * nElements[1] + BasisFunctionType::nDofsPerNode()) * dofsPerRow;
  
  int elementZ = int(elementNo / (nElements[0] * nElements[1]));
  int elementY = int((elementNo % (nElements[0] * nElements[1])) / nElements[0]);
  int elementX = elementNo % nElements[0];
  int localZ = int(dofIndex / MathUtility::sqr(BasisFunctionType::nDofsPerBasis()));
  int localY = int((dofIndex % MathUtility::sqr(BasisFunctionType::nDofsPerBasis())) / BasisFunctionType::nDofsPerBasis());
  int localX = dofIndex % BasisFunctionType::nDofsPerBasis();
  
  return dofsPerPlane * (elementZ * averageNDofsPerElement1D + localZ)
    + dofsPerRow * (elementY * averageNDofsPerElement1D + localY) 
    + averageNDofsPerElement1D * elementX + localX;
}

//! get all dofs of a specific node for 3D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<int> dofGlobalNos) const
{
  for (int i=0; i<BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos.push_back(BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i);
  }
}

// element-local nodeIndex to global nodeNo for 1D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0 1
  // averageNDofsPerElement: 
  // 1         2            2
  // nDofsPerBasis:
  // 2         3            4
  // nDofsPerNode:
  // 1         1            2
  // nNodesPerElement:
  // 2         3            2
  return BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNNodesPerElement() * elementNo + nodeIndex;
}

// element-local nodeIndex to global nodeNo for 2D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      2 3
  // 0 1       0 1 2      0 1
  // nNodesPerElement:
  // 4         9          4
  
  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given
  const std::array<element_no_t,2> nElements{this->nElements(0), this->nElements(1)};
  
  int averageNNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement();
  int nodesPerRow = (averageNNodesPerElement1D * nElements[0] + 1);
  int elementX = int(elementNo % nElements[0]);
  int elementY = int(elementNo / nElements[0]);
  int localX = nodeIndex % nNodesPerElement1D;
  int localY = int(nodeIndex / nNodesPerElement1D);
  
  return nodesPerRow * (elementY * averageNNodesPerElement1D + localY) 
    + averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex to global nodeNo for 3D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given
  const std::array<element_no_t,3> nElements{this->nElements(0), this->nElements(1), this->nElements(2)};
  
  int averageNNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement();
  int nodesPerRow = (averageNNodesPerElement1D * nElements[0] + 1);
  int nodesPerPlane = (averageNNodesPerElement1D * nElements[1] + 1) * nodesPerRow;
  
  int elementZ = int(elementNo / (nElements[0] * nElements[1]));
  int elementY = int((elementNo % (nElements[0] * nElements[1])) / nElements[0]);
  int elementX = elementNo % nElements[0];
  int localZ = int(nodeIndex / MathUtility::sqr(nNodesPerElement1D));
  int localY = int((nodeIndex % MathUtility::sqr(nNodesPerElement1D)) / nNodesPerElement1D);
  int localX = nodeIndex % nNodesPerElement1D;
  
  return nodesPerPlane * (elementZ * averageNNodesPerElement1D + localZ)
    + nodesPerRow * (elementY * averageNNodesPerElement1D + localY) 
    + averageNNodesPerElement1D * elementX + localX;
}


template<int D,typename BasisFunctionType>
element_no_t BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
nElements() const
{
  return this->nElements_;
}

template<int D,typename BasisFunctionType>
BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
BasisOnMeshDofs(PyObject *settings, bool noGeometryField) :
  BasisOnMeshJacobian<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::BasisOnMeshJacobian(settings),
  noGeometryField_(noGeometryField)
{
  LOG(TRACE) << "BasisOnMeshDofs constructor";
 
  if (PythonUtility::containsKey(settings, "exelem"))
  {
    std::string filenameExelem = PythonUtility::getOptionString(settings, "exelem", "input.exelem");
    std::string filenameExnode = PythonUtility::getOptionString(settings, "exnode", "input.exnode");
    
    // ------------------------------------------------------------------------
    // read in exelem file
    this->parseExelemFile(filenameExelem);
    
    // ------------------------------------------------------------------------
    // read in exnode file
    this->parseExnodeFile(filenameExnode);
    
    // remap names of field variables if specified in config
    this->remapFieldVariables(settings);
    
    
    // eliminate scale factors
    //this->eliminateScaleFactors();
  }
  else if (PythonUtility::containsKey(settings, "nodePositions"))
  {
    this->parseFromSettings(settings);
  }
  else
  {
    LOG(FATAL) << "Could not create UnstructuredDeformable node positions. " 
      << "Either specify \"exelem\" and \"exnode\" or \"nodePositions\". ";
  }
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
parseFromSettings(PyObject *settings)
{
  LOG(TRACE) << "parseFromSettings";
 
  // example input in settings:
  //  "nodePositions": [[0,0,0], [1,0], [2,0,0], [0,1], [1,1], [2,1], [0,2], [1,2], [2,2], ...],
  
  // parse node positions
  std::vector<Vec3> nodePositions;
  
  // get the first node position from the list
  PyObject *pyNodePositions = PythonUtility::getOptionListBegin<PyObject *>(settings, "nodePositions");

  // loop over other entries of list
  for (;
      !PythonUtility::getOptionListEnd(settings, "nodePositions");
      PythonUtility::getOptionListNext<PyObject *>(settings, "nodePositions", pyNodePositions))
  {
    Vec3 nodePosition = PythonUtility::convertFromPython<Vec3>(pyNodePositions);
    LOG(DEBUG) << "node position " << nodePosition;
    nodePositions.push_back(nodePosition);
  }
  
  // parse elements
  struct Element
  {
    struct ElementNode
    {
      node_no_t nodeGlobalNo;
      unsigned int versionNo;
    };
    std::vector<ElementNode> nodes;
  };
  std::vector<Element> elements;
  
  // example input in settings:
  //  "elements": [[[0,0], [1,0], [2,1], [3,0]], [next element]]   # each node is [node no, version-at-that-node no] or just node-no then it assumes version no 0
  
  // get the first node position from the list
  PyObject *pyElement = PythonUtility::getOptionListBegin<PyObject *>(settings, "elements");

  // loop over other entries of list
  for (;
      !PythonUtility::getOptionListEnd(settings, "elements");
      PythonUtility::getOptionListNext<PyObject *>(settings, "elements", pyElement))
  {
    // get the python list that makes up the element, e.g. [[0,0], [1,0], [2,1], [3,0]]
    typedef std::array<PyObject *,this->nNodesPerElement()> PyElementNodes;
    PyElementNodes pyElementNodes = PythonUtility::convertFromPython<PyObject *,this->nNodesPerElement()>(pyElement,PyLong_FromLong(this->nNodesPerElement()));
   
    Element currentElement;
    currentElement.nodes.resize(this->nNodesPerElement());
    // loop over nodes of that element
    for (int nodeIndex = 0; nodeIndex < this->nNodesPerElement(); nodeIndex++)
    {
       // extract the node positions, e.g. [1,0] (global node no., version no.) or just 1 (only global node no., version no. defaults to 0)
       std::array<int,2> elementNode = PythonUtility::convertFromPython<int,2>(pyElementNodes[nodeIndex], 0);
       
       LOG(DEBUG) << "   elementNode " << elementNode;
       
       currentElement.nodes[nodeIndex].nodeGlobalNo = elementNode[0];
       currentElement.nodes[nodeIndex].versionNo = elementNode[1];
    }
    elements.push_back(currentElement);
  }
  
  this->nElements_ = elements.size();
  
  LOG(DEBUG) << nodePositions.size() << " nodes, " << elements.size() << " elements";
  
  // initialize elementToNodeMapping_
  this->elementToNodeMapping_ = std::make_shared<FieldVariable::ElementToNodeMapping>();
  this->elementToNodeMapping_->setNumberElements(elements.size());
  
  // loop over elements
  for (element_no_t elementGlobalNo = 0; elementGlobalNo < elements.size(); elementGlobalNo++)
  {
    // allocate vector in elementToNodeMapping for nodes of that element
    this->elementToNodeMapping_->getElement(elementGlobalNo).nodeGlobalNo.resize(this->nNodesPerElement());
    
    // loop over nodes of that element
    for (int nodeIndex = 0; nodeIndex < this->nNodesPerElement(); nodeIndex++)
    {
      // assign previously extracted global node number to element in ElementToNodeMapping
      this->elementToNodeMapping_->getElement(elementGlobalNo).nodeGlobalNo[nodeIndex] = elements[elementGlobalNo].nodes[nodeIndex].nodeGlobalNo;
      
      // check if given global node number is valid
      if (this->elementToNodeMapping_->getElement(elementGlobalNo).nodeGlobalNo[nodeIndex] >= nodePositions.size()) 
      {
        LOG(FATAL) << "Element " << elementGlobalNo << " contains node global no. "
          << this->elementToNodeMapping_->getElement(elementGlobalNo).nodeGlobalNo[nodeIndex] << " which is >= the number of nodes (" 
          << nodePositions.size() << ")";
      }
    }
  }
  
  // set exfileRepresentation
  std::shared_ptr<FieldVariable::ExfileRepresentation> exfileRepresentation = std::make_shared<FieldVariable::ExfileRepresentation>();
  exfileRepresentation->setNumberElements(elements.size());
  
  // loop over elements
  for (element_no_t elementGlobalNo = 0; elementGlobalNo < elements.size(); elementGlobalNo++)
  {
    std::shared_ptr<FieldVariable::ExfileElementRepresentation> &exfileElementRepresentation 
     = exfileRepresentation->getExfileElementRepresentation(elementGlobalNo);
    
    // construct exfile element representation for the current element 
    exfileElementRepresentation = std::make_shared<FieldVariable::ExfileElementRepresentation>();
    exfileElementRepresentation->setNumberNodes(this->nNodesPerElement());
    
    // loop over nodes of current element 
    for (int nodeIndex = 0; nodeIndex < this->nNodesPerElement(); nodeIndex++)
    {
      FieldVariable::ExfileElementRepresentation::Node &node = exfileElementRepresentation->getNode(nodeIndex);
      /*
      struct Node
      {
        std::vector<int> valueIndices;        ///< the indices of the dof values of this node in the exnode file node values (sub-)block (for the particular field variable/component) of the node. If there are not multiple versions, this is simply 0,1,...,ndofs-1. If there are e.g. 2 versions and 8 dofs per node, this can be 0,1,...,7 if the elements uses the 1st version, or 8,...,15 if the element uses the second version. Note, that the real index of the dofs inside the values block may be different when this is not the first component of the block.
        std::vector<int> scaleFactorIndices;   ///< the indices of all scale factor entries for this node in the exelem element scale factors block. Thus this is kind of a node to element block mapping. 
      };*/
      
      // set valueIndices
      int versionNo = elements[elementGlobalNo].nodes[nodeIndex].versionNo;
      node.valueIndices.resize(this->nDofsPerNode());
      std::iota(node.valueIndices.begin(), node.valueIndices.end(), versionNo*this->nDofsPerNode());
    }
  }
  
  // setup elementToDof mapping, this also creates nodeToDofMapping
  std::shared_ptr<FieldVariable::ElementToDofMapping> elementToDofMapping = std::make_shared<FieldVariable::ElementToDofMapping>();
  
  elementToDofMapping->setNumberElements(elements.size());
  std::shared_ptr<FieldVariable::NodeToDofMapping> nodeToDofMapping 
    = elementToDofMapping->setup(exfileRepresentation, this->elementToNodeMapping_, this->nDofsPerNode());
  
  LOG(DEBUG) << "nodeToDofMapping: " << *nodeToDofMapping;
    
  this->nDofs_ = elementToDofMapping->nDofs();
    
  // setup geometric field variable
  std::pair<std::string, std::shared_ptr<FieldVariableType>> geometryField(
    "geometry", 
    std::make_shared<FieldVariableType>());
  this->fieldVariable_.insert(geometryField);
  
  this->fieldVariable_["geometry"]->initializeFromMappings("geometry", true, exfileRepresentation, elementToDofMapping, 
                                                           this->elementToNodeMapping_, nodeToDofMapping, {"x","y","z"});
  
  // unify exfile representation variables in field variables such that there is only one shared for all components of a field variable
  this->fieldVariable_["geometry"]->unifyMappings(this->elementToNodeMapping_, this->nDofsPerNode());
  
  // set values from nodePositions
  
  // loop over nodes
  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < nodePositions.size(); nodeGlobalNo++)
  { 
    int nVersions = nodeToDofMapping->nVersions(nodeGlobalNo);
   
    LOG(DEBUG) << "node " << nodeGlobalNo << ", nVersions: " << nVersions;
    
    // loop over components
    int componentNo = 0;
    auto &component = this->fieldVariable_["geometry"]->component();
    for (auto iter = component.begin(); iter != component.end(); iter++, componentNo++)
    {
      // get value of nodePositions vector that matches current component
      double positionValue = nodePositions[nodeGlobalNo][componentNo];
      
      // create vector containing all dofs of the current node
      std::vector<double> nodeValues(this->nDofsPerNode()*nVersions, 0.0);
      
      // set first dof of every version for the particular component (this leaves derivative dofs of Hermite at 0)
      for (int versionNo = 0; versionNo < nVersions; versionNo++)
      {
        LOG(DEBUG) << "    set version " << versionNo << " at " << this->nDofsPerNode()*versionNo << " to " << positionValue;
        nodeValues[this->nDofsPerNode()*versionNo] = positionValue;
      }

      LOG(DEBUG) << "   component " << iter->first << ", positionValue: " << positionValue << " nodeValues: " <<nodeValues;
      
      // set nodal dof values at node
      iter->second.setNodeValues(nodeGlobalNo, nodeValues.begin());
    }
  }
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
initialize()
{
  if (this->noGeometryField_)
    return;
  
  std::shared_ptr<BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> ptr = this->shared_from_this();
  
  assert(ptr != nullptr);
  
  std::shared_ptr<BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> self = std::static_pointer_cast<BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(ptr);
  
  assert(self != nullptr);
  assert(fieldVariable_.find("geometry") != fieldVariable_.end());
  
  this->fieldVariable_.at("geometry")->setMesh(self);
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
remapFieldVariables(PyObject *settings)
{
  // remap name of field variables to different names if specified
  if (PythonUtility::containsKey(settings, "remap"))
  {
    std::string keyString = "remap";
    std::pair<std::string, std::string> dictItem 
      = PythonUtility::getOptionDictBegin<std::string, std::string>(settings, keyString);
    
    for (; !PythonUtility::getOptionDictEnd(settings, keyString); 
        PythonUtility::getOptionDictNext<std::string, std::string>(settings, keyString, dictItem))
    {
      std::string key = dictItem.first;
      std::string value = dictItem.second;
          
      if (this->fieldVariable_.find(key) != this->fieldVariable_.end())
      {
        std::shared_ptr<FieldVariableType> fieldVariable = this->fieldVariable_[key];
        this->fieldVariable_.erase(key);
        this->fieldVariable_[value] = fieldVariable;
      }
      else
      {
        LOG(WARNING) << "Remap of field variable from \"" << key << "\" to \"" << value << "\" failed: no such entry.";
      }
    }
  }
  
  // if there is no field with name "geometry"
  if (this->fieldVariable_.find("geometry") == this->fieldVariable_.end())
  { 
    bool geometryFieldFound = false;
    // search for a geometry field 
    for(auto &fieldVariableEntry : this->fieldVariable_)
    {
      if (fieldVariableEntry.second->isGeometryField())
      {
        LOG(WARNING) << "Remap geometry field variable from \"" << fieldVariableEntry.first << "\" to \"geometry\".";
      
        std::shared_ptr<FieldVariableType> fieldVariable = fieldVariableEntry.second;
        this->fieldVariable_.erase(fieldVariableEntry.first);
        this->fieldVariable_["geometry"] = fieldVariable;
        
        geometryFieldFound = true;
        break;
      }
    }
    
    // output field variables
    for(auto &fieldVariableEntry : this->fieldVariable_)
    {
      VLOG(1) << *fieldVariableEntry.second;
    }
    
    if (!geometryFieldFound)
    {
      LOG(FATAL) << "The specified Exfiles contain no geometry field. The field must be named \"geometry\" or have the type \"coordinates\".";
    }
  }
}

template<int D,typename BasisFunctionType>
int BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  if (this->fieldVariable_.find("geometry") == this->fieldVariable_.end())
    LOG(FATAL) << "Mesh contains no field variable \"geometry\". Use remap to create one!";
  
  return this->fieldVariable_.at("geometry")->getDofNo(elementNo, dofIndex);
}

template<int D,typename BasisFunctionType>
int BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  return this->elementToNodeMapping_->getElement(elementNo).nodeGlobalNo[nodeIndex];
}

//! get all dofs of a specific node, unstructured mesh
template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<int> dofGlobalNos) const
{
  if (this->fieldVariable_.find("geometry") == this->fieldVariable_.end())
    LOG(FATAL) << "Mesh contains no field variable \"geometry\". Use remap to create one!";
  
  std::vector<int> &nodeDofs = dofGlobalNos.insert(dofGlobalNos.end(), 
    this->fieldVariable_.at("geometry")->nodeToDofMapping()->getNodeDofs(nodeGlobalNo));
  dofGlobalNos.reserve(dofGlobalNos.size() + nodeDofs.size());
  
  for(int dof : nodeDofs)
  {
    dofGlobalNos.push_back(dof);
  }
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
eliminateScaleFactors()
{
  // loop over field variables
  for (auto fieldVariable : this->fieldVariable_)
  {
    fieldVariable.second->eliminateScaleFactors();
  }
}

template<int D,typename BasisFunctionType>
dof_no_t BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
nDofs() const
{
  return nDofs_;
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
addNonGeometryFieldVariables(std::vector<std::shared_ptr<FieldVariableType>> &fieldVariables)
{
  // loop over field variables
  for (auto fieldVariable : this->fieldVariable_)
  {
    if (!fieldVariable.second->isGeometryField())
      fieldVariables.push_back(fieldVariable.second);
  }
}
  
};  // namespace

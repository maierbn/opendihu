#include "basis_on_mesh/03_basis_on_mesh_dofs.h"

#include <cmath>
#include <array>
#include <string>
#include <map>
#include <cassert>

#include "easylogging++.h"
#include "utility/string_utility.h"
#include "utility/math_utility.h"

#include "field_variable/unstructured/exfile_representation.h"
#include "field_variable/unstructured/element_to_dof_mapping.h"

namespace BasisOnMesh
{
 
using namespace StringUtility;

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
 
  if (PythonUtility::hasKey(settings, "exelem"))
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
    
    // eliminate scale factors (not yet tested)
    //this->eliminateScaleFactors();
  }
  else if (PythonUtility::hasKey(settings, "nodePositions"))
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
    VLOG(1) << "node position " << nodePosition;
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
       
       VLOG(1) << "   elementNode " << elementNode;
       
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
  
  VLOG(1) << "nodeToDofMapping: " << *nodeToDofMapping;
    
  this->nDofs_ = elementToDofMapping->nDofs();
    
  // create and setup geometric field variable
  this->geometryField_ = std::make_shared<FieldVariable::FieldVariable<BasisOnMeshType,3>>();
  
  this->geometryField_->initializeFromMappings("geometry", true, exfileRepresentation, elementToDofMapping, 
                                                           this->elementToNodeMapping_, nodeToDofMapping, {"x","y","z"});
  
  // unify exfile representation variables in field variables such that there is only one shared for all components of a field variable
  this->geometryField_->unifyMappings(this->elementToNodeMapping_, this->nDofsPerNode());
  
  // set values from nodePositions
  
  // loop over nodes
  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < nodePositions.size(); nodeGlobalNo++)
  { 
    int nVersions = nodeToDofMapping->nVersions(nodeGlobalNo);
   
    VLOG(1) << "node " << nodeGlobalNo << ", nVersions: " << nVersions;
    
    // loop over components
    int componentNo = 0;
    std::array<::FieldVariable::Component<BasisOnMeshType>,3> &component = this->geometryField_->component();
    for (auto iter = component.begin(); iter != component.end(); iter++, componentNo++)
    {
      // get value of nodePositions vector that matches current component
      double positionValue = nodePositions[nodeGlobalNo][componentNo];
      
      // create vector containing all dofs of the current node
      std::vector<double> nodeValues(this->nDofsPerNode()*nVersions, 0.0);
      
      // set first dof of every version for the particular component (this leaves derivative dofs of Hermite at 0)
      for (int versionNo = 0; versionNo < nVersions; versionNo++)
      {
        VLOG(1) << "    set version " << versionNo << " at " << this->nDofsPerNode()*versionNo << " to " << positionValue;
        nodeValues[this->nDofsPerNode()*versionNo] = positionValue;
      }

      VLOG(1) << "   component no. " << componentNo << ", positionValue: " << positionValue << " nodeValues: " <<nodeValues;
      
      // set nodal dof values at node
      iter->setNodeValues(nodeGlobalNo, nodeValues.begin());
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
  assert(geometryField_ != nullptr);
  
  this->geometryField_->setMesh(self);
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
remapFieldVariables(PyObject *settings)
{
  // remap name of field variables to different names if specified
  if (PythonUtility::hasKey(settings, "remap"))
  {
    std::string keyString = "remap";
    std::pair<std::string, std::string> dictItem 
      = PythonUtility::getOptionDictBegin<std::string, std::string>(settings, keyString);
    
    for (; !PythonUtility::getOptionDictEnd(settings, keyString); 
        PythonUtility::getOptionDictNext<std::string, std::string>(settings, keyString, dictItem))
    {
      // remap request field variable from key to value
      std::string key = dictItem.first;
      std::string value = dictItem.second;
          
      if (this->fieldVariable_.find(key) != this->fieldVariable_.end())
      {
        std::shared_ptr<FieldVariableBaseType> fieldVariable = this->fieldVariable_[key];
        this->fieldVariable_.erase(key);
        this->fieldVariable_[value] = fieldVariable;
      }
      else
      {
        LOG(WARNING) << "Remap of field variable from \"" << key << "\" to \"" << value << "\" failed: no such entry.";
      }
    }
  }
  
  // if there is no geometry field stored yet, extract from named variables
  if (!this->geometryField_)
  {
    if (this->fieldVariable_.find("geometry") != this->fieldVariable_.end())
    {
      this->geometryField_ = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,3>>(
         this->fieldVariable_.at("geometry"));
      this->fieldVariable_.erase("geometry");
    }
  }
  
  // if there is still no geometry field stored
  if (!this->geometryField_)
  { 
    bool geometryFieldFound = false;
    
    // search for a geometry field that is named differently
    for(auto &fieldVariableEntry : this->fieldVariable_)
    {
      if (fieldVariableEntry.second->isGeometryField())
      {
        LOG(WARNING) << "Remap geometry field variable from \"" << fieldVariableEntry.first << "\" to \"geometry\".";
      
        std::shared_ptr<FieldVariableBaseType> fieldVariable = fieldVariableEntry.second;
        this->fieldVariable_.erase(fieldVariableEntry.first);
        this->geometryField_ = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,3>>(fieldVariable);
        
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
dof_no_t BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  if (!this->geometryField_)
    LOG(FATAL) << "Mesh contains no field variable \"geometry\". Use remap to create one!";
  
  return this->geometryField_->getDofNo(elementNo, dofIndex);
}

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  return this->elementToNodeMapping_->getElement(elementNo).nodeGlobalNo[nodeIndex];
}

//! get all dofs of a specific node, unstructured mesh
template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  if (!this->geometryField_)
    LOG(FATAL) << "Mesh contains no field variable \"geometry\". Use remap to create one!";
  
  std::vector<dof_no_t> &nodeDofs = this->geometryField_->nodeToDofMapping()->getNodeDofs(nodeGlobalNo);
  
  dofGlobalNos.reserve(dofGlobalNos.size() + nodeDofs.size());
  
  for(dof_no_t dof : nodeDofs)
  {
    dofGlobalNos.push_back(dof);
  }
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
eliminateScaleFactors()
{
  // this method was not tested yet!
  // loop over field variables
  for (auto fieldVariable : this->fieldVariable_)
    fieldVariable.second->eliminateScaleFactors();
  
  if (geometryField_)
    geometryField_->eliminateScaleFactors();
}

template<int D,typename BasisFunctionType>
dof_no_t BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
nDofs() const
{
  return nDofs_;
}

template<int D,typename BasisFunctionType>
std::shared_ptr<::FieldVariable::FieldVariableBase<BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> 
BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
fieldVariable(std::string name)
{
  if (fieldVariable_.find(name) != fieldVariable_.end())
    return fieldVariable_.at(name);
  else
    return nullptr;
}
  
};  // namespace

#include "function_space/04_function_space_data_unstructured.h"

#include "easylogging++.h"
#include "utility/string_utility.h"

#include "basis_function/basis_function.h"
#include "field_variable/factory.h"

#include <iostream>
#include <fstream>

namespace FunctionSpace
{

template<int D,typename BasisFunctionType>
void FunctionSpaceDataUnstructured<D,BasisFunctionType>::
parseFromSettings(PyObject *settings)
{
  LOG(TRACE) << "parseFromSettings";

  // example input in settings:
  //  "nodePositions": [[0,0,0], [1,0], [2,0,0], [0,1], [1,1], [2,1], [0,2], [1,2], [2,2], ...],

  // parse node positions
  std::vector<Vec3> nodePositions;
  bool listWarningIssued = false;    // if the warning about lists was already shown

  // get the first node position from the list
  PyObject *pyNodePositions = PythonUtility::getOptionListBegin<PyObject *>(settings, "nodePositions");

  // loop over other entries of list
  for (;
      !PythonUtility::getOptionListEnd(settings, "nodePositions");
      PythonUtility::getOptionListNext<PyObject *>(settings, "nodePositions", pyNodePositions))
  {
    if (!PythonUtility::isTypeList(pyNodePositions) && !listWarningIssued)
    {
      listWarningIssued = true;
      LOG(WARNING) << "\"nodePositions\" is not a list of lists.";
    }
   
    Vec3 nodePosition = PythonUtility::convertFromPython<Vec3>::get(pyNodePositions);
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
    PyElementNodes pyElementNodes = PythonUtility::convertFromPython<PyElementNodes>::
      get(pyElement,PyLong_FromLong(this->nNodesPerElement()));

    Element currentElement;
    currentElement.nodes.resize(this->nNodesPerElement());
    // loop over nodes of that element
    for (int nodeIndex = 0; nodeIndex < this->nNodesPerElement(); nodeIndex++)
    {
       // extract the node positions, e.g. [1,0] (global node no., version no.) or just 1 (only global node no., version no. defaults to 0)
       std::array<int,2> elementNode = PythonUtility::convertFromPython<std::array<int,2>>::get(pyElementNodes[nodeIndex], 0);

       VLOG(1) << "   elementNode " << elementNode;

       currentElement.nodes[nodeIndex].nodeGlobalNo = elementNode[0];
       currentElement.nodes[nodeIndex].versionNo = elementNode[1];
    }
    elements.push_back(currentElement);
  }

  this->nElements_ = elements.size();

  LOG(DEBUG) << nodePositions.size() << " node positions, " << elements.size() << " elements";

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

  this->nDofs_ = elementToDofMapping->nDofsLocal();

  // create and setup geometry field variable
  this->geometryField_ = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>();
  
  
  // retrieve "this" pointer and convert to downwards pointer of most derived class "FunctionSpace"
  std::shared_ptr<FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> thisFunctionSpace
    = std::static_pointer_cast<FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(this->shared_from_this());

  assert(thisFunctionSpace != nullptr);
  assert(geometryField_ != nullptr);
  
  // set this mesh to geometry field
  geometryField_->setFunctionSpace(thisFunctionSpace);

  this->geometryField_->initializeFromMappings("geometry", true, exfileRepresentation, elementToDofMapping,
                                                           this->elementToNodeMapping_, nodeToDofMapping, {"x","y","z"});
  //this does not yet call initializeValuesVector(), which needs to be done after the meshPartition of the mesh is set
  
  // unify exfile representation variables in field variables such that there is only one shared for all components of a field variable
  this->geometryField_->unifyMappings(this->elementToNodeMapping_, this->nDofsPerNode());

  // create meshPartition and redistribute elements if necessary, this needs information about mesh size
  FunctionSpacePartition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::initialize();
  
  // create the values vector
  this->geometryField_->initializeValuesVector();
  
  // set values from nodePositions

  // loop over nodes
  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < nodePositions.size(); nodeGlobalNo++)
  {
    int nVersions = nodeToDofMapping->nVersions(nodeGlobalNo);  // this fails

    VLOG(1) << "node " << nodeGlobalNo << ", nVersions: " << nVersions;

    // loop over components
    int componentNo = 0;
    std::array<::FieldVariable::Component<FunctionSpaceType,3>,3> &component = this->geometryField_->component();
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

};  // namespace

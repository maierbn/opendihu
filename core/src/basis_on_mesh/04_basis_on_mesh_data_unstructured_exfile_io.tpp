#include "basis_on_mesh/04_basis_on_mesh_data_unstructured.h"

#include "easylogging++.h"
#include "utility/string_utility.h"

#include "basis_function/basis_function.h"
#include "field_variable/factory.h"

#include <iostream>
#include <fstream>

namespace BasisOnMesh
{

template<int D,typename BasisFunctionType>
void BasisOnMeshDataUnstructured<D,BasisFunctionType>::
parseExelemFile(std::string exelemFilename)
{
  std::ifstream file_exelem(exelemFilename.c_str(), std::ios::in | std::ios::binary);
  if (!file_exelem.is_open())
  {
    LOG(WARNING) << "Could not open exelem file \"" << exelemFilename << "\" for reading.";
  }

  VLOG(1) << "parseExelemFile";

  // first pass of file: find out number of elements in file
  this->nElements_ = 0;
  std::string line;
  while(!file_exelem.eof())
  {
    getline(file_exelem, line);
    if(file_exelem.eof())
      break;

    if (line.find("Element:") != std::string::npos)
    {
      element_no_t elementGlobalNo = getNumberAfterString(line, "Element:");
      this->nElements_ = std::max(this->nElements_, elementGlobalNo);
    }
  }

  if (this->elementToNodeMapping_ == nullptr)
    this->elementToNodeMapping_ = std::make_shared<FieldVariable::ElementToNodeMapping>();
  this->elementToNodeMapping_->setNumberElements(this->nElements_);

  VLOG(1) << "nElements: " <<this->nElements_;

  // reset file stream
  file_exelem.clear();
  file_exelem.seekg(0, std::ios::beg);

  // second pass of file: read in dofs for each element
  int fieldNo = 0;
  //int nFields;

  enum
  {
    nothing,
    elementsFollow,     // after "Shape. Dimension=3" was found
    fieldsFollow,       // after "#Fields"
  }
  lineType = nothing;

  std::stringstream nextFieldNo;
  std::string fieldName;
  std::string fieldContent;
  std::string elementContent;

  std::vector<int> valueIndices, scaleFactorIndices;

  // loop over lines of file
  while(!file_exelem.eof())
  {
    // parse line
    getline(file_exelem, line);
    if(file_exelem.eof())
      break;

    // check if line contains "Shape."
    if (line.find("Shape.") != std::string::npos && line.find("Dimension=") != std::string::npos)
    {
      int dimension = getNumberAfterString(line, "Dimension=");

      VLOG(1) << "dimension=" <<dimension;
      if (dimension == D)
      {
        lineType = elementsFollow;
        VLOG(1) << "elementsFollow";
      }
      else if(dimension > D)
      {
        LOG(WARNING) << "Exfile contains " << dimension << "D data, UnstructuredDeformable mesh is " << D << "D.";
      }
      continue;
    }

    // check if line contains "#Fields="
    if(lineType == elementsFollow && line.find("#Fields=") != std::string::npos)
    {
      // parse number of following fields, not used
      //nFields = getNumberAfterString(line, "#Fields=");

      // prepare string that starts the next field description block
      fieldNo = 0;
      nextFieldNo.str("");
      nextFieldNo << fieldNo+1 << ")";
      lineType = fieldsFollow;
      VLOG(1) << "fieldsFollow";

      // finish last element
      if (!elementContent.empty())
      {
        // parse whole element block into elementToNode mapping
        this->elementToNodeMapping_->parseElementFromExelemFile(elementContent);
        for(auto &fieldVariable : this->fieldVariable_)
        {
          fieldVariable.second->parseElementFromExelemFile(elementContent);
        }
        elementContent = "";
      }

      continue;
    }

    // if in field description block
    if (lineType == fieldsFollow)
    {
      // if the next field description block begins or a block with elements begins (i.e. the current field description block is finished)
      if(line.find(nextFieldNo.str()) != std::string::npos || line.find("Element:") != std::string::npos)
      {

        VLOG(1) << "next field description block begins";

        // finish previous fieldVariable component
        if (fieldNo != 0)
        {
          VLOG(1) << "finish previous fieldVariable component of name [" << fieldName << "]";

          // extract number of components from fieldContent
          if (fieldContent.find("#Components=") == std::string::npos)
          {
            LOG(ERROR) << "Could not parse number of components for field variable \"" << fieldName << "\".";
          }
          const int nComponents = StringUtility::getNumberAfterString(fieldContent, "#Components=");

          // if the field variable with this name does not exist already, create new field variable object
          if (this->fieldVariable_.find(fieldName) == this->fieldVariable_.end())
          {
            std::pair<std::string, std::shared_ptr<FieldVariableBaseType>> newEntry
            (
              fieldName,
              FieldVariable::Factory<BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::makeShared(nComponents)
            );
            this->fieldVariable_.insert(newEntry);
          }
          this->fieldVariable_[fieldName]->setNumberElements(this->nElements_);

          VLOG(1) << "new field variable [" << fieldName << "]";

          // parse whole field description block inside component object
          this->fieldVariable_[fieldName]->parseHeaderFromExelemFile(fieldContent);
          fieldContent = "";
        }

        // get name of current fieldVariable, if there is one starting here
        if (line.find(nextFieldNo.str()) != std::string::npos)
        {
          std::string name = line;
          extractUntil(name, nextFieldNo.str());
          fieldName = extractUntil(name, ",");
          trim(fieldName);
        }

        // increase number of next field variable and prepare string to search for for next field variable description block
        fieldNo++;
        nextFieldNo.str("");
        nextFieldNo << fieldNo+1 << ")";
      }

      // if a block with elements begins at current line
      if (line.find("Element:") != std::string::npos)
      {
        lineType = elementsFollow;
        elementContent = "";
      }
      else
      {
        // collect lines of current field description block
        fieldContent += line+"\n";
      }
    }

    // if the current line is part of an elements block
    if (lineType == elementsFollow)
    {
      // if a new element begins, i.e. the previous element is finished
      if (line.find("Element:") != std::string::npos)
      {
        VLOG(1) << "new element";
        // if there were lines collected for the previous element block
        if (!elementContent.empty())
        {
          // parse whole element block into elementToNode mapping
          this->elementToNodeMapping_->parseElementFromExelemFile(elementContent);
          for(auto &fieldVariable : this->fieldVariable_)
          {
            fieldVariable.second->parseElementFromExelemFile(elementContent);
          }
          elementContent = "";
        }
      }

      // collect lines for current element block
      elementContent += line+"\n";
    }
  }

  // if there were lines collected for the previous element block
  if (!elementContent.empty())
  {
    // parse whole element block into elmeentToNode mapping
    this->elementToNodeMapping_->parseElementFromExelemFile(elementContent);
    for(auto &fieldVariable : this->fieldVariable_)
    {
      fieldVariable.second->parseElementFromExelemFile(elementContent);
    }
  }

  VLOG(1) << "parsing done";
  for(auto &fieldVariable : this->fieldVariable_)
  {
    VLOG(1)<<*fieldVariable.second;
  }

  VLOG(1) << "unifyMappings - 1";
  // unify exfile representation variables in field variables such that there is only one shared for all components of a field variable
  for(auto &fieldVariable : this->fieldVariable_)
  {
    fieldVariable.second->unifyMappings(this->elementToNodeMapping_, this->nDofsPerNode());
  }

  VLOG(1) << "unifyMappings - 2";
  // remove duplicate exfile representation and element to dof mapping entries that are shared among different field variables
  for(auto &fieldVariable : this->fieldVariable_)
  {
    for(auto &fieldVariable2 : this->fieldVariable_)
    {
      if (fieldVariable.first != fieldVariable2.first)
      {
        VLOG(1) << " unifyMappings(fieldVariable) " << fieldVariable.first << "," << fieldVariable2.first;
        fieldVariable.second->unifyMappings(fieldVariable2.second);
      }
    }
  }

  // get number of dofs
  this->nDofs_ = 0;

  if (!this->fieldVariable_.empty())
  {
    std::shared_ptr<FieldVariable::ElementToDofMapping> elementToDofMapping = this->fieldVariable_.begin()->second->elementToDofMapping();
    this->nDofs_ = elementToDofMapping->nLocalDofs();
  }
  VLOG(1) << "nDofs: " << this->nDofs_;

  //! check if

  VLOG(1) << "initialize values vectors";
  // allocate common values vector for each field variable
  for(auto &fieldVariable : this->fieldVariable_)
  {
    fieldVariable.second->initializeValuesVector();
  }
  VLOG(1) << "parseExelemFile done";
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDataUnstructured<D,BasisFunctionType>::
parseExnodeFile(std::string exnodeFilename)
{
  std::ifstream file_exnode(exnodeFilename.c_str(), std::ios::in | std::ios::binary);
  if (!file_exnode.is_open())
  {
    LOG(WARNING) << "Could not open exnode file \"" << exnodeFilename << "\" for reading.";
  }

  VLOG(1) << "parseExnodeFile";
  // read in file content
  std::string content( (std::istreambuf_iterator<char>(file_exnode) ),
                       (std::istreambuf_iterator<char>()    ) );

  for(auto &fieldVariable : this->fieldVariable_)
  {
    // set all values to 0.0
    fieldVariable.second->setValues(0.0);

    // parse file and values0
    fieldVariable.second->parseFromExnodeFile(content);
  }
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDataUnstructured<D,BasisFunctionType>::
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

  this->nDofs_ = elementToDofMapping->nLocalDofs();

  // create and setup geometry field variable
  this->geometryField_ = std::make_shared<FieldVariable::FieldVariable<BasisOnMeshType,3>>();

  this->geometryField_->initializeFromMappings("geometry", true, exfileRepresentation, elementToDofMapping,
                                                           this->elementToNodeMapping_, nodeToDofMapping, {"x","y","z"});

  // unify exfile representation variables in field variables such that there is only one shared for all components of a field variable
  this->geometryField_->unifyMappings(this->elementToNodeMapping_, this->nDofsPerNode());

  // set values from nodePositions

  // loop over nodes
  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < nodePositions.size(); nodeGlobalNo++)
  {
    int nVersions = nodeToDofMapping->nVersions(nodeGlobalNo);  // this fails

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
void BasisOnMeshDataUnstructured<D,BasisFunctionType>::
initializeGeometryField()
{
  // this gets called from initialize in 06_basis_on_mesh_dofs_nodes
  // it has to determine the number of elements or everything that is needed to create a meshPartition 
  // (this is done in initialize of 03_basis_on_mesh_partition_unstructured.tpp)
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDataUnstructured<D,BasisFunctionType>::
setGeometryFieldValues()
{
  // this gets called from initialize in 06_basis_on_mesh_dofs_nodes
  // after the geometry field got their mesh and meshPartition
 
}

};  // namespace

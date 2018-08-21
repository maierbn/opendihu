#include "field_variable/unstructured/02_field_variable_data_unstructured_deformable.h"

#include <sstream>
#include <map>
#include <memory>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "utility/string_utility.h"

namespace FieldVariable
{

using namespace StringUtility;

// contructor as data copy with a different name (component names are the same)
template<int D, typename BasisFunctionType, int nComponents>
FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
FieldVariableData(FieldVariable<BasisOnMeshType,nComponents> &rhs, std::string name) :
  FieldVariableComponents<BasisOnMeshType,nComponents>::FieldVariableComponents()
{
  // create new distributed petsc vec as copy of rhs values vector
  this->values_ = std::make_shared<PartitionedPetscVec<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>>(
     rhs.partitionedPetscVec(), name);

  // initialize everything from other field variable
  initializeFromFieldVariable(rhs, name, rhs.componentNames());
}

// constructor with mesh, name and components
template<int D, typename BasisFunctionType, int nComponents>
FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
FieldVariableData(std::shared_ptr<BasisOnMeshType> mesh, std::string name, std::vector<std::string> componentNames) :
  FieldVariableComponents<BasisOnMeshType,nComponents>::FieldVariableComponents()
{
  assert(false); // not implemented
  // this is needed for mixed formulation, to implement this and set as the lower order field variable I need the higher order field variable as parameter
}

template<int D, typename BasisFunctionType, int nComponents>
FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
FieldVariableData() :
  FieldVariableComponents<BasisOnMeshType,nComponents>::FieldVariableComponents()
{
}

template<int D, typename BasisFunctionType, int nComponents>
FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
~FieldVariableData()
{
  if (this->values_)
  {
    //PetscErrorCode ierr = VecDestroy(&*this->values_); CHKERRV(ierr);
  }
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
setMesh(std::shared_ptr<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> mesh)
{
  this->mesh_ = mesh;
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
parseHeaderFromExelemFile(std::string content)
{
  VLOG(2) << "parseHeaderFromExelemFile(" << content << ")";
  bool newComponentStarts = false;
  int nComponentsParsed = 0;
  int componentNo = 0;
  int nLinesOfNodesFollow = 0;

  std::string componentName = "";
  std::string componentContent = "";

  // loop over content line-wise
  unsigned int pos = 0;
  while(pos < content.size())
  {
    // extract next line
    unsigned int posNewline = content.find("\n",pos);
    std::string line = content.substr(pos, posNewline-pos);
    if (posNewline == std::string::npos)
      pos = content.size();
    else
      pos = posNewline+1;

    // if new component description in header starts
    if (line.find(")") != std::string::npos)
    {
      trim(line);
      std::string no = extractUntil(line,")");
      trim(no);
      exfileNo_ = atoi(no.c_str());

      this->name_ = extractUntil(line, ",");
      trim(this->name_);
      std::string type = extractUntil(line, ",");
      this->isGeometryField_ = type.find("coordinate") != std::string::npos;
      nComponentsParsed = getNumberAfterString(line, "#Components=");

      if (nComponentsParsed != nComponents)
      {
        LOG(DEBUG) << "line=[" << line << "]";
        LOG(ERROR) << "Number of components is not parsed correctly. nComponentsParsed=" << nComponentsParsed << ", nComponents=" << nComponents;
      }

      newComponentStarts = true;
      VLOG(2) << "nComponents: " << nComponentsParsed;

      componentNo = 0;
      continue;
    }

    // new block of component definition starts
    if (newComponentStarts)
    {
      // finish previous component
      if (componentNo != 0)
      {

        VLOG(2) << "finish previous component, parse header";
        component_[componentNo-1].initialize(nullptr, componentNo-1, nElements_);
        component_[componentNo-1].parseHeaderFromExelemFile(componentContent);
        this->componentNames_[componentNo-1] = componentName;
        componentContent = "";
      }
      componentNo++;
      componentContent += line+"\n";
      componentName = extractUntil(line, ".");
      trim(componentName);
      newComponentStarts = false;
      continue;
    }

    // number of nodes
    if (line.find("#Nodes=") != std::string::npos)
    {
      componentContent += line+"\n";
      int nNodes = getNumberAfterString(line, "#Nodes=");
      VLOG(2) << "nNodes: " << nNodes;
      nLinesOfNodesFollow = nNodes*3;
      continue;
    }

    // further component block
    if (nLinesOfNodesFollow > 0)
    {
      componentContent += line+"\n";
      nLinesOfNodesFollow--;
      if (nLinesOfNodesFollow == 0)
      {
        newComponentStarts = true;
      }
    }
  }

  // finish previous component
  if (componentNo != 0)
  {
    VLOG(2) << "finish previous component, parse header";

    component_[componentNo-1].initialize(nullptr, componentNo-1, nElements_);
    component_[componentNo-1].parseHeaderFromExelemFile(componentContent);
    this->componentNames_[componentNo-1] = componentName;
  }
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
parseElementFromExelemFile(std::string content)
{
  for (auto &component : component_)
  {
    component.parseElementFromExelemFile(content);
  }
}
/*
template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
setNumberElements(element_no_t nElements)
{
  this->nElements_ = nElements;
}

template<int D, typename BasisFunctionType, int nComponents>
element_no_t FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
nElementsLocal() const
{
  return this->nElements_;
}*/

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
unifyMappings(std::shared_ptr<ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode)
{
  VLOG(1) << "unifyMappings 1 with " << component_.size() << " components";

  // unify exfileRepresentation object
  // check if exfileRepresentation is equal for all components
  int componentNo = 0;
  if (componentNo < nComponents)
  {
    VLOG(1) << " component " << this->componentNames_[componentNo];

    Component<BasisOnMeshType,nComponents> &firstComponent = component_[componentNo];
    std::string firstComponentName = this->componentNames_[componentNo];

    // extract exfile representation from first component
    this->exfileRepresentation_ = firstComponent.exfileRepresentation();
    assert(this->exfileRepresentation_);

    // loop over further components after the first component
    componentNo++;
    for (; componentNo < nComponents; componentNo++)
    {
      VLOG(1) << " component " << this->componentNames_[componentNo];
      Component<BasisOnMeshType,nComponents> &component = component_[componentNo];

      assert(component.exfileRepresentation());

      if (*this->exfileRepresentation_ == *component.exfileRepresentation())
      {
        VLOG(1) << " set new exfile representation";
        component.setExfileRepresentation(this->exfileRepresentation_);
      }
      else
      {
        LOG(ERROR) << "The components \"" << firstComponentName << "\" and \"" << this->componentNames_[componentNo]
          << "\" of field variable \"" << this->name_ << "\" have different exfile representations.";
      }
    }
  }

  createElementToDofMapping(elementToNodeMapping, nDofsPerNode);
}

template<int D, typename BasisFunctionType, int nComponents>
template<int nComponents2>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
unifyMappings(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents2> &fieldVariable)
{
  VLOG(1) << "unifyMappings 2";
  // loop over own components
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    Component<BasisOnMeshType,nComponents> &component = component_[componentNo];
    VLOG(1) << "first: " << this->componentNames_[componentNo];

    // loop over components of the other fieldVariable
    for (int componentNo2 = 0; componentNo2 < nComponents2; componentNo2++)
    //for (typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter2 = fieldVariable.component_.begin(); iter2 != fieldVariable.component_.end(); iter2++)
    {
      Component<BasisOnMeshType,nComponents2> &component2 = fieldVariable.component_[componentNo2];
      VLOG(1) << "second: " << fieldVariable.componentNames_[componentNo2];

      // assert that pointers are not null
      assert (component.exfileRepresentation());
      assert (component2.exfileRepresentation());
      assert (component.elementToDofMapping());
      assert (component2.elementToDofMapping());

      if (*component.exfileRepresentation() == *component2.exfileRepresentation())
      {
        VLOG(1) << "set exfile rep for " << fieldVariable.componentNames_[componentNo2] << " to be the same as for " << this->componentNames_[componentNo];
        component2.setExfileRepresentation(component.exfileRepresentation());
      }
      else VLOG(1) << "exfileRepresentation is different";

      if (*component.elementToDofMapping() == *component2.elementToDofMapping())
      {
        VLOG(1) << "set elementToDof and nodeToDofMapping for " << fieldVariable.componentNames_[componentNo2] << " to be the same as for " << this->componentNames_[componentNo];
        component2.setDofMappings(component.elementToDofMapping(), component.nodeToDofMapping());
      }
      else VLOG(1) << "elementToDofMapping is different";
    }
  }
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
unifyMappings(std::shared_ptr<FieldVariableBase<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> fieldVariable2)
{
  VLOG(1) << "unifyMappings 2";
  // loop over own components
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    Component<BasisOnMeshType,nComponents> &component = component_[componentNo];
    VLOG(1) << "first: " << this->componentNames_[componentNo];

    // loop over components of the other fieldVariable
    for (int componentNo2 = 0; componentNo2 < fieldVariable2->getNComponents(); componentNo2++)
    {
      // depending on number of components of the fieldVariable2, get different type of Component<...,nComponents2>
      if (fieldVariable2->getNComponents() == 1)
      {
        const int nComponents2 = 1;
        std::shared_ptr<FieldVariableComponents<BasisOnMeshType,nComponents2>> fieldVariableComponents2
          = std::static_pointer_cast<FieldVariableComponents<BasisOnMeshType,nComponents2>>(fieldVariable2);
        std::shared_ptr<Component<BasisOnMeshType,nComponents2>> component2 = fieldVariableComponents2->component(componentNo2);
      
        if (!component2)
        {
          continue;
        }
        VLOG(1) << "second: " << fieldVariable2->componentName(componentNo2);

        // assert that pointers are not null
        assert (component.exfileRepresentation());
        assert (component2->exfileRepresentation());
        assert (component.elementToDofMapping());
        assert (component2->elementToDofMapping());

        if (*component.exfileRepresentation() == *(component2->exfileRepresentation()))
        {
          VLOG(1) << "set exfile rep for " << fieldVariable2->componentName(componentNo2) << " to be the same as for " << this->componentNames_[componentNo];
          component2->setExfileRepresentation(component.exfileRepresentation());
        }
        else VLOG(1) << "exfileRepresentation is different";

        if (*component.elementToDofMapping() == *(component2->elementToDofMapping()))
        {
          VLOG(1) << "set elementToDof and nodeToDofMapping for " << fieldVariable2->componentName(componentNo2)
            << " to be the same as for " << this->componentNames_[componentNo];
          component2->setDofMappings(component.elementToDofMapping(), component.nodeToDofMapping());
        }
        else VLOG(1) << "elementToDofMapping is different";
      }
      else if (fieldVariable2->getNComponents() == 2)
      {
        const int nComponents2 = 2;
        std::shared_ptr<FieldVariableComponents<BasisOnMeshType,nComponents2>> fieldVariableComponents2
          = std::static_pointer_cast<FieldVariableComponents<BasisOnMeshType,nComponents2>>(fieldVariable2);
        std::shared_ptr<Component<BasisOnMeshType,nComponents2>> component2 = fieldVariableComponents2->component(componentNo2);
      
        if (!component2)
        {
          continue;
        }
        VLOG(1) << "second: " << fieldVariable2->componentName(componentNo2);

        // assert that pointers are not null
        assert (component.exfileRepresentation());
        assert (component2->exfileRepresentation());
        assert (component.elementToDofMapping());
        assert (component2->elementToDofMapping());

        if (*component.exfileRepresentation() == *(component2->exfileRepresentation()))
        {
          VLOG(1) << "set exfile rep for " << fieldVariable2->componentName(componentNo2) << " to be the same as for " << this->componentNames_[componentNo];
          component2->setExfileRepresentation(component.exfileRepresentation());
        }
        else VLOG(1) << "exfileRepresentation is different";

        if (*component.elementToDofMapping() == *(component2->elementToDofMapping()))
        {
          VLOG(1) << "set elementToDof and nodeToDofMapping for " << fieldVariable2->componentName(componentNo2)
            << " to be the same as for " << this->componentNames_[componentNo];
          component2->setDofMappings(component.elementToDofMapping(), component.nodeToDofMapping());
        }
        else VLOG(1) << "elementToDofMapping is different";
      }
      else if (fieldVariable2->getNComponents() == 3)
      {
        const int nComponents2 = 3;
        std::shared_ptr<FieldVariableComponents<BasisOnMeshType,nComponents2>> fieldVariableComponents2
          = std::static_pointer_cast<FieldVariableComponents<BasisOnMeshType,nComponents2>>(fieldVariable2);
        std::shared_ptr<Component<BasisOnMeshType,nComponents2>> component2 = fieldVariableComponents2->component(componentNo2);
      
        if (!component2)
        {
          continue;
        }
        VLOG(1) << "second: " << fieldVariable2->componentName(componentNo2);

        // assert that pointers are not null
        assert (component.exfileRepresentation());
        assert (component2->exfileRepresentation());
        assert (component.elementToDofMapping());
        assert (component2->elementToDofMapping());

        if (*component.exfileRepresentation() == *(component2->exfileRepresentation()))
        {
          VLOG(1) << "set exfile rep for " << fieldVariable2->componentName(componentNo2) << " to be the same as for " << this->componentNames_[componentNo];
          component2->setExfileRepresentation(component.exfileRepresentation());
        }
        else VLOG(1) << "exfileRepresentation is different";

        if (*component.elementToDofMapping() == *(component2->elementToDofMapping()))
        {
          VLOG(1) << "set elementToDof and nodeToDofMapping for " << fieldVariable2->componentName(componentNo2)
            << " to be the same as for " << this->componentNames_[componentNo];
          component2->setDofMappings(component.elementToDofMapping(), component.nodeToDofMapping());
        }
        else VLOG(1) << "elementToDofMapping is different";
      }
      else if (fieldVariable2->getNComponents() == 4)
      {
        const int nComponents2 = 4;
        std::shared_ptr<FieldVariableComponents<BasisOnMeshType,nComponents2>> fieldVariableComponents2
          = std::static_pointer_cast<FieldVariableComponents<BasisOnMeshType,nComponents2>>(fieldVariable2);
        std::shared_ptr<Component<BasisOnMeshType,nComponents2>> component2 = fieldVariableComponents2->component(componentNo2);
      
        if (!component2)
        {
          continue;
        }
        VLOG(1) << "second: " << fieldVariable2->componentName(componentNo2);

        // assert that pointers are not null
        assert (component.exfileRepresentation());
        assert (component2->exfileRepresentation());
        assert (component.elementToDofMapping());
        assert (component2->elementToDofMapping());

        if (*component.exfileRepresentation() == *(component2->exfileRepresentation()))
        {
          VLOG(1) << "set exfile rep for " << fieldVariable2->componentName(componentNo2) << " to be the same as for " << this->componentNames_[componentNo];
          component2->setExfileRepresentation(component.exfileRepresentation());
        }
        else VLOG(1) << "exfileRepresentation is different";

        if (*component.elementToDofMapping() == *(component2->elementToDofMapping()))
        {
          VLOG(1) << "set elementToDof and nodeToDofMapping for " << fieldVariable2->componentName(componentNo2)
            << " to be the same as for " << this->componentNames_[componentNo];
          component2->setDofMappings(component.elementToDofMapping(), component.nodeToDofMapping());
        }
        else VLOG(1) << "elementToDofMapping is different";
      }
      else 
      {
        LOG(WARNING) << "unifyMappings not implemented for field variables with > 4 components (this is not an error, it only means that more memory is needed).";
      }
    }
  }
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
createElementToDofMapping(std::shared_ptr<ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode)
{
  if (!this->elementToDofMapping_)
    this->elementToDofMapping_ = std::make_shared<ElementToDofMapping>();

  // set number of elements in element to dof mapping
  this->elementToDofMapping_->setNumberElements(this->nElements_);

  // create the element to dof mapping from exfileRepresentation and element to node mapping
  this->nodeToDofMapping_ = this->elementToDofMapping_->setup(this->exfileRepresentation_, elementToNodeMapping, nDofsPerNode);

  // set the element to dof and node to dof mapping at each component
  for (auto &component : component_)
  {
    component.setDofMappings(this->elementToDofMapping_, this->nodeToDofMapping_);
  }
  elementToNodeMapping_ = elementToNodeMapping;
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
parseFromExnodeFile(std::string content)
{
  //int nFields;
  int fieldNo = 0;
  std::stringstream nextFieldNo;
  std::map<std::string, int> componentDataBlockStartIndex;
  std::vector<double> blockValues;
  int nodeNo = 0;
  int previousNodeNo = 0;

  enum {
    nothing,
    fieldsFollow,
    componentsFollow,
    valuesFollow,
  } lineType = nothing;

  // loop over file content line-wise
  unsigned int pos = 0;
  while(pos < content.size())
  {
    // extract next line
    unsigned int posNewline = content.find("\n",pos);
    std::string line = content.substr(pos, posNewline-pos);
    if (posNewline == std::string::npos)
      pos = content.size();
    else
      pos = posNewline+1;

    VLOG(2) << "[" << StringUtility::replace(line, "\r", "") << "]";

    bool finishPreviousNodeDataBlock = false;

    // line with number of fields, e.g. "#Fields=2", this indicates starting field description block
    if (line.find("#Fields=") != std::string::npos)
    {
      // finish previous values block
      if (lineType == valuesFollow)
      {
        // finish previous node
        previousNodeNo = nodeNo;
        finishPreviousNodeDataBlock = true;
      }

      // parse number of fields in current field description block, not used
      //nFields = getNumberAfterString(line, "#Fields=");

      // prepare string to look for for next field description
      fieldNo = 0;
      nextFieldNo.str("");
      nextFieldNo << fieldNo+1 << ") ";

      VLOG(2) << "search for [" << nextFieldNo.str() << "]";
      lineType = fieldsFollow;
    }

    // line with field, e.g. "1) coordinates, coordinate, rectangular cartesian, #Component<BasisOnMeshType>s=3"
    else if ((lineType == fieldsFollow || lineType == componentsFollow)
      && line.find(nextFieldNo.str()) != std::string::npos)
    {
      // extract name of field
      extractUntil(line, ")");
      std::string name = extractUntil(line, ",");
      trim(name);

      VLOG(2) << "field \"" << name << "\".";

      // consider it if it matches the own name
      if (name == this->name_)
      {
        lineType = componentsFollow;
        VLOG(2) << "matches own name, parse components";
      }
      else
      {
        lineType = fieldsFollow;
      }

      // prepare string to search for for the next field
      fieldNo++;
      nextFieldNo.str("");
      nextFieldNo << fieldNo+1 << ") ";
      VLOG(2) << "search for [" << nextFieldNo.str() << "]";
      continue;
    }

    // line with beginning of node, e.g. "Node: 73"
    if (line.find("Node:") != std::string::npos)
    {
      previousNodeNo = nodeNo;
      nodeNo = getNumberAfterString(line, "Node:")-1;
      VLOG(2) << "  node " << nodeNo;

      if (lineType == valuesFollow)
      {
        // finish previous node
        finishPreviousNodeDataBlock = true;
      }
      else
      {
        lineType = valuesFollow;
        continue;
      }
    }

    // line with component, e.g. "x.  Value index=1, #Derivatives=7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions=2"
    if (lineType == componentsFollow)
    {
      std::string componentName = extractUntil(line, ".");
      trim(componentName);

      // search for component name
      bool componentNameFound = false;
      for (auto component : this->componentNames_)
      {
        if (component == componentName)
        {
          componentNameFound = true;
          break;
        }
      }

      if (!componentNameFound)
      {
        LOG(WARNING) << "Component \"" << componentName << "\" of field variable \"" << this->name_ << "\" in exnode file is not present in exelem file.";
      }
      else
      {
        componentDataBlockStartIndex[componentName] = getNumberAfterString(line, "index=") - 1;
        VLOG(2) << "component name " << componentName << " index=" << componentDataBlockStartIndex[componentName];
      }
    }

    // if a node values block ended prior to the current line, store the collected values for the component
    if (finishPreviousNodeDataBlock)
    {
      VLOG(2) << "finishPreviousNodeDataBlock";
      // loop over components of current field variable
      for (auto &pair : componentDataBlockStartIndex)
      {
        std::string componentName = pair.first;
        int subBlockStartIndex = pair.second;
        int componentNo = this->findComponent(componentName);

        VLOG(2) << "  componentName=[" << componentName<<"], subBlockStartIndex=" << subBlockStartIndex
          << " n values: " << blockValues.size() << " previousNodeNo=" << previousNodeNo;

        // store data at component
        component_[componentNo].setNodeValuesFromBlock(previousNodeNo, blockValues.begin()+subBlockStartIndex);
      }
      // clear temporary data vector
      blockValues.clear();
      continue;
    }

    // values
    if (lineType == valuesFollow)
    {
      VLOG(2) << "valuesFollow";
      trim(line);
      while(!line.empty())
      {
        blockValues.push_back(atof(line.c_str()));
        // proceed to next number
        if (line.find(" ") != std::string::npos)
        {
          extractUntil(line, " ");
          trim(line);
        }
        else break;
      }
      VLOG(2) << "extract values block " << blockValues;
    }
  }   // while


  if (lineType == valuesFollow)
  {
    // loop over components of current field variable
    for (auto &pair : componentDataBlockStartIndex)
    {
      // get data block index of component
      std::string componentName = pair.first;
      int subBlockStartIndex = pair.second;
      int componentNo = this->findComponent(componentName);

      // store data at component
      component_[componentNo].setNodeValuesFromBlock(nodeNo, blockValues.begin()+subBlockStartIndex);
    }
    // clear temporary data vector
    blockValues.clear();
  }

  // finalize Petsc vectors
  this->finishVectorManipulation();
}

template<int D, typename BasisFunctionType, int nComponents>
dof_no_t FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  return this->component_.begin()->getDofNo(elementNo, dofIndex);
}

template<int D, typename BasisFunctionType, int nComponents>
std::shared_ptr<Component<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>> FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
component(int componentNo)
{
  assert(componentNo >= 0);
  assert(componentNo < nComponents);
  return std::make_shared<Component<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>>(this->component_[componentNo]);
}

template<int D, typename BasisFunctionType, int nComponents>
Component<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> &FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
component(std::string componentName)
{
  int componentNo = this->componentNames_.find(componentName);
  return this->component_[componentNo];
}

template<int D, typename BasisFunctionType, int nComponents>
std::array<Component<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>,nComponents> &FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
component()
{
  return this->component_;
}

template<int D, typename BasisFunctionType, int nComponents>
std::shared_ptr<ExfileRepresentation> FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
exfileRepresentation() const
{
  return this->exfileRepresentation_;
}

template<int D, typename BasisFunctionType, int nComponents>
std::shared_ptr<ElementToDofMapping> FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
elementToDofMapping() const
{
  return this->elementToDofMapping_;
}

template<int D, typename BasisFunctionType, int nComponents>
std::shared_ptr<ElementToNodeMapping> FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
elementToNodeMapping() const
{
  return this->elementToNodeMapping_;
}

template<int D, typename BasisFunctionType, int nComponents>
std::shared_ptr<NodeToDofMapping> FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
nodeToDofMapping() const
{
  return this->nodeToDofMapping_;
}

template<int D, typename BasisFunctionType, int nComponents>
Vec &FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
valuesLocal()
{
  return this->values_->valuesLocal();
}

template<int D, typename BasisFunctionType, int nComponents>
Vec &FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
valuesGlobal()
{
  return this->values_->valuesGlobal();
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
initializeValuesVector()
{
  // initialize the PETSc vector that contains all the value entries
  // create PartitionedPetscVector, the size of the vector is stored in meshPartition
  if (this->values_ == nullptr)
  {
    this->values_ = std::make_shared<PartitionedPetscVec<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>>(
      this->mesh_->meshPartition(), this->name_
    );
  }

  // set vector for all components
  for (auto &component : this->component_)
  {
    component.setValuesVector(this->values_);
  }
}/*

template<int D, typename BasisFunctionType, int nComponents>
std::size_t FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
nEntries() const
{
  return this->nEntries_;
}

template<int D, typename BasisFunctionType, int nComponents>
dof_no_t FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
nDofsLocalWithGhosts() const
{
  return this->nEntries_ / nComponents;
}

template<int D, typename BasisFunctionType, int nComponents>
node_no_t FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
nNodesLocalWithGhosts() const
{
  // because parallelism is not implemented for unstructured meshes, there is no difference to with or without ghosts
  return nodeToDofMapping_->nLocalNodes();
}
*/
template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
initializeComponents(std::vector<std::string> &componentNames, std::string exfileBasisRepresentation)
{
  // insert components
  int componentIndex = 0;

  assert(componentNames.size() == nComponents);
  std::copy(componentNames.begin(), componentNames.end(), this->componentNames_.begin());

  // create a new values vector for the new field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    Component<BasisOnMeshType,nComponents> &component = this->component_[componentNo];
    component.initialize(this->values_, componentIndex++, this->nElements_);   // note: this->values_ may be nullptr but is updated by initializeValuesVector
    component.setName(this->componentNames_[componentNo], exfileBasisRepresentation);
    component.setDofMappings(this->elementToDofMapping_, this->nodeToDofMapping_);
    component.setExfileRepresentation(this->exfileRepresentation_);
  }

  // set values_ and nEntries_
  initializeValuesVector();
}

template<int D, typename BasisFunctionType, int nComponents>
template<typename FieldVariableType>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
initializeFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames)
{
  // set member variables
  this->name_ = name;
  this->exfileNo_ = 0;
  this->nEntries_ = 0;
  this->isGeometryField_ = false;
  this->nElements_ = fieldVariable.nElementsLocal();
  this->exfileRepresentation_ = fieldVariable.exfileRepresentation();
  this->elementToDofMapping_ = fieldVariable.elementToDofMapping();
  this->elementToNodeMapping_ = fieldVariable.elementToNodeMapping();
  this->nodeToDofMapping_ = fieldVariable.nodeToDofMapping();
  this->mesh_ = fieldVariable.mesh();

  std::string exfileBasisRepresentation = fieldVariable.component().begin()->exfileBasisFunctionSpecification();

  initializeComponents(componentNames, exfileBasisRepresentation);
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
initializeFromMappings(std::string name, bool isGeometryField,
                       std::shared_ptr<ExfileRepresentation> exfileRepresentation,
                       std::shared_ptr<ElementToDofMapping> elementToDofMapping,
                       std::shared_ptr<ElementToNodeMapping> elementToNodeMapping,
                       std::shared_ptr<NodeToDofMapping> nodeToDofMapping,
                       std::vector<std::string> componentNames)
{
  this->name_ = name;
  this->exfileNo_ = 0;
  this->nEntries_ = 0;
  this->isGeometryField_ = isGeometryField;
  this->nElements_ = elementToDofMapping->nElementsLocal();

  // remove duplicate exfile representations
  exfileRepresentation->unifyExfileElementRepresentations();

  this->exfileRepresentation_ = exfileRepresentation;
  this->elementToDofMapping_ = elementToDofMapping;
  this->elementToNodeMapping_ = elementToNodeMapping;
  this->nodeToDofMapping_ = nodeToDofMapping;

  // this->mesh still needs to be set by setMesh

  std::string exfileBasisRepresentation = BasisFunction::getBasisRepresentationString<D,BasisFunctionType>();

  initializeComponents(componentNames, exfileBasisRepresentation);


  LOG(DEBUG) << "FieldVariable nDofs: " << this->nDofsLocal();
}

template<int D, typename BasisFunctionType, int nComponents>
int FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
getNumberScaleFactors(element_no_t elementGlobalNo) const
{

  //! return the node numbers and scale factors of the element
  return elementToNodeMapping_->getElement(elementGlobalNo).scaleFactors.size();
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
outputHeaderExelem(std::ostream &stream, element_no_t currentElementGlobalNo, int fieldVariableNo)
{
  // set no of field variable in ex file if it was specified
  if (fieldVariableNo != -1)
    exfileNo_ = fieldVariableNo+1;

  // output first line of header
  stream << " " << exfileNo_ << ") " << this->name_ << ", " << (this->isGeometryField_? "coordinate" : "field")
    << ", rectangular cartesian, #Components=" << nComponents << std::endl;

  // output headers of components
  for (auto &component : component_)
  {
    component.outputHeaderExelem(stream, currentElementGlobalNo);
  }
  /*
   1) Geometry, coordinate, rectangular cartesian, #Components=2
     x.   l.Lagrange*l.Lagrange, no modify, standard node based.
     #Nodes= 4
     1.  #Values=1
      Value indices:     1
      Scale factor indices:    1
     2.  #Values=1
      Value indices:     1
      Scale factor indices:    2
   */
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
outputHeaderExnode(std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo)
{
  // set no of field variable in ex file if it was specified
  if (fieldVariableNo != -1)
    exfileNo_ = fieldVariableNo+1;

  // output first line of header
  stream << " " << exfileNo_ << ") " << this->name_ << ", " << (this->isGeometryField_? "coordinate" : "field")
    << ", rectangular cartesian, #Components=" << nComponents << std::endl;

  // output headers of components
  for (auto &component : component_)
  {
    component.outputHeaderExnode(stream, currentNodeGlobalNo, valueIndex);
  }
}

template<int D, typename BasisFunctionType, int nComponents>
bool FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
{
  // loop over components
  for (auto &component : component_)
  {
    if(!component.exfileRepresentation()->haveSameExfileRepresentation(element1, element2))
    {
      LOG(DEBUG) << "  component " << component << " has different exfileRepr for elements " << element1 << " and " << element2;
      return false;
    }
  }
  return true;
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
eliminateScaleFactors()
{
  // loop over elements
  for (element_no_t elementGlobalNo = 0; elementGlobalNo < nElements_; elementGlobalNo++)
  {
    // loop over components
    for (auto &component : component_)
    {
      ElementToNodeMapping::Element &element = elementToNodeMapping_->getElement(elementGlobalNo);
      // loop over element dofs
      for (unsigned int nodeIdx = 0; nodeIdx < element.nodeGlobalNo.size(); nodeIdx++)
      {
        node_no_t nodeGlobalNo = element.nodeGlobalNo[nodeIdx];
        double scaleFactor = element.scaleFactors[nodeIdx];

        std::vector<int> &nodeDofs = nodeToDofMapping_->getNodeDofs(nodeGlobalNo);
        std::vector<double> &scaleFactors = nodeToDofMapping_->getNodeScaleFactors(nodeGlobalNo);
        std::vector<double> nodeValues(nodeDofs.size());

        // loop over node dofs
        for (unsigned int dofIdx = 0; dofIdx < nodeDofs.size(); dofIdx++)
        {
          // get dof value
          double value = component.getValue(nodeDofs[dofIdx]);

          double scaleFactor2 = scaleFactors[dofIdx];
          assert(scaleFactor - scaleFactor2 < 1e-16);

          LOG(DEBUG) << "scaleFactor: " << scaleFactor << "," << scaleFactor2;

          // multiply value with scale factor
          value *= scaleFactor;

          nodeValues[dofIdx] = value;
        }

        // set updated values
        component.setValuesForNode(nodeGlobalNo, nodeValues.begin());
      }
    }
  }
}

template<int D, typename BasisFunctionType, int nComponents>
void FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
output(std::ostream &stream) const
{
  stream << "\"" << this->name_ << "\", nElements: " << nElements_
    << ", isGeometryField: " << std::boolalpha << this->isGeometryField_ << std::endl
    << "  components:" << std::endl;
  for (auto &componentName : this->componentNames_)
  {
    stream << componentName;
  }
  stream << "  exfileRepresentation: " << std::endl;
  if (exfileRepresentation_ == nullptr)
    stream << "null" << std::endl;
  else
    stream << *exfileRepresentation_ << std::endl;
  /*
    << "  elementToDofMapping: " << std::endl
    << *elementToDofMapping_ << std::endl
    << "  elementToNodeMapping: " << std::endl
    << *elementToNodeMapping_ << std::endl
    << "  nodeToDofMapping: " << std::endl
    << *nodeToDofMapping_ << std::endl;
  */
}

template<int D, typename BasisFunctionType, int nComponents>
std::shared_ptr<PartitionedPetscVec<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>>
FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
partitionedPetscVec()
{
  return values_; 
}

};  // namespace

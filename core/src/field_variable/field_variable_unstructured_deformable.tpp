#include "field_variable/field_variable_unstructured_deformable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <map>
#include <fstream>
#include <iomanip>

namespace FieldVariable
{
  
using namespace StringUtility;
  
template<int D, typename BasisFunctionType>
FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
~FieldVariable()
{
  if (this->values_)
  {
    //PetscErrorCode ierr = VecDestroy(&*this->values_); CHKERRV(ierr);
  }
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
parseHeaderFromExelemFile(std::string content)
{
  VLOG(2) << "parseHeaderFromExelemFile(" << content << ")";
  bool newComponentStarts = false;
  int nComponents = 0;
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
      isGeometryField_ = type.find("coordinate") != std::string::npos;
      nComponents = getNumberAfterString(line, "#Components=");
      newComponentStarts = true;
      VLOG(2) << "nComponents: " << nComponents;
      
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
       
        component_[componentName].initialize(nullptr, nComponents, componentNo-1, nElements_);
        component_[componentName].parseHeaderFromExelemFile(componentContent);
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
   
    component_[componentName].initialize(nullptr, nComponents, componentNo-1, nElements_);
    component_[componentName].parseHeaderFromExelemFile(componentContent);
  }
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
parseElementFromExelemFile(std::string content)
{
  for (typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter = component_.begin(); iter != component_.end(); iter++)
  {
    Component<BasisOnMeshType> &component = iter->second;
    component.parseElementFromExelemFile(content);
  }
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
setNumberElements(element_no_t nElements)
{
  this->nElements_ = nElements;
}

template<int D, typename BasisFunctionType>
element_no_t FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nElements() const
{
  return this->nElements_;
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
unifyMappings(std::shared_ptr<ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode)
{
 
  VLOG(1) << "unifyMappings 1 with " << component_.size() << " components";
 
  // unify exfileRepresentation object
  // check if exfileRepresentation is equal for all components
  typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter = component_.begin();
  if (iter != component_.end())
  { 
    VLOG(1) << " component " << iter->first;
   
    Component<BasisOnMeshType> &firstComponent = iter->second;
    std::string firstComponentKey = iter->first;
    
    // extract exfile representation from first component
    this->exfileRepresentation_ = firstComponent.exfileRepresentation();
    assert(this->exfileRepresentation_);
    
    // loop over further component after the first component
    iter++;
    for (; iter != component_.end(); iter++)
    {
      VLOG(1) << " component " << iter->first;
      Component<BasisOnMeshType> &component = iter->second;
    
      assert(component.exfileRepresentation());
      
      if (*this->exfileRepresentation_ == *component.exfileRepresentation())
      {
        VLOG(1) << " set new exfile representation";
        component.setExfileRepresentation(this->exfileRepresentation_);
      }
      else 
      {
        LOG(ERROR) << "The components \"" << firstComponentKey << "\" and \"" << iter->first 
          << "\" of field variable \"" << this->name_ << "\" have different exfile representations.";
      }
    }
  }
  
  createElementToDofMapping(elementToNodeMapping, nDofsPerNode);
}
  
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
unifyMappings(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> &fieldVariable)
{
  VLOG(1) << "unifyMappings 2";
  // loop over own components
  for (typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter = component_.begin(); iter != component_.end(); iter++)
  {
    Component<BasisOnMeshType> &component = iter->second;
    VLOG(1) << "first: " << iter->first;
    
    // loop over components of the other fieldVariable
    for (typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter2 = fieldVariable.component_.begin(); iter2 != fieldVariable.component_.end(); iter2++)
    {
      Component<BasisOnMeshType> &component2 = iter2->second;
      VLOG(1) << "second: " << iter2->first;
      
      // assert that pointers are not null
      assert (component.exfileRepresentation());
      assert (component2.exfileRepresentation());
      assert (component.elementToDofMapping());
      assert (component2.elementToDofMapping());
      
      if (*component.exfileRepresentation() == *component2.exfileRepresentation())
      {
        VLOG(1) << "set exfile rep for " << iter2->first << " to be the same as for " << iter->first;
        component2.setExfileRepresentation(component.exfileRepresentation());
      }
      else VLOG(1) << "exfileRepresentation is different";
      
      if (*component.elementToDofMapping() == *component2.elementToDofMapping())
      {
        VLOG(1) << "set elementToDof and nodeToDofMapping for " << iter2->first << " to be the same as for " << iter->first;
        component2.setDofMappings(component.elementToDofMapping(), component.nodeToDofMapping());
      }
      else VLOG(1) << "elementToDofMapping is different";
    }
  }
  
}
  
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
createElementToDofMapping(std::shared_ptr<ElementToNodeMapping> elementToNodeMapping, const int nDofsPerNode)
{
  if (!this->elementToDofMapping_)
    this->elementToDofMapping_ = std::make_shared<ElementToDofMapping>();
 
  // set number of elements in element to dof mapping
  this->elementToDofMapping_->setNumberElements(this->nElements_);
  
  // create the element to dof mapping from exfileRepresentation and element to node mapping
  this->nodeToDofMapping_ = this->elementToDofMapping_->setup(this->exfileRepresentation_, elementToNodeMapping, nDofsPerNode);
  
  // set the element to dof and node to dof mapping at each component
  for (typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter = component_.begin(); iter != component_.end(); iter++)
  {
    Component<BasisOnMeshType> &component = iter->second;
    component.setDofMappings(this->elementToDofMapping_, this->nodeToDofMapping_);
  }
  elementToNodeMapping_ = elementToNodeMapping;
}
  
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
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
      
      if (component_.find(componentName) == component_.end())
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
        
        VLOG(2) << "  componentName=[" << componentName<<"], subBlockStartIndex=" << subBlockStartIndex 
          << " n values: " << blockValues.size() << " previousNodeNo=" << previousNodeNo;
        
        // store data at component
        component_[componentName].setNodeValuesFromBlock(previousNodeNo, blockValues.begin()+subBlockStartIndex);
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
      
      // store data at component
      component_[componentName].setNodeValuesFromBlock(nodeNo, blockValues.begin()+subBlockStartIndex);
    }
    // clear temporary data vector
    blockValues.clear();
  }
}

template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  return this->component_.begin()->second.getDofNo(elementNo, dofIndex);
}

template<int D, typename BasisFunctionType>
Component<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> &FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
component(std::string key)
{
  return this->component_[key]; 
}
  
template<int D, typename BasisFunctionType>
std::map<std::string, Component<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>> &FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
component()
{
  return this->component_; 
}

template<int D, typename BasisFunctionType>
std::shared_ptr<ExfileRepresentation> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
exfileRepresentation() const
{
  return this->exfileRepresentation_; 
}

template<int D, typename BasisFunctionType>
std::shared_ptr<ElementToDofMapping> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
elementToDofMapping() const
{
  return this->elementToDofMapping_; 
}

template<int D, typename BasisFunctionType>
std::shared_ptr<ElementToNodeMapping> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
elementToNodeMapping() const
{
  return this->elementToNodeMapping_; 
}

template<int D, typename BasisFunctionType>
std::shared_ptr<NodeToDofMapping> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nodeToDofMapping() const
{
  return this->nodeToDofMapping_;
}

template<int D, typename BasisFunctionType>
Vec &FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
values()
{
  return *this->values_; 
}

//! for a specific component, get values from their global dof no.s
template<int D, typename BasisFunctionType>
template<int N>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getValues(std::string component, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values)
{
  this->component_[component].getValues(dofGlobalNo, values);
}
/*
template<int D, typename BasisFunctionType>
template<int N, int nComponents>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
{
  std::array<double,nComponents> resultVector;
  
  // transform global dof no.s to vector indices of first component
  for (int valueIndex = 0; valueIndex < N; valueIndex++)
  {
    int valuesVectorIndex = dofGlobalNo[valueIndex]*nComponents;
    
    // create indices vector with values {0,1,2,...,nComponents-1}
    std::array<int,nComponents> indices;
    for(int i=0; i<nComponents; i++)
      indices[i] = valuesVectorIndex + i;
    
    // get values and assign them to result values vector
    VecGetValues(*this->values_, nComponents, indices.data(), resultVector.data());
    values[valueIndex] = resultVector;
  }
}*/

/*
//! for a specific component, get the values corresponding to all element-local dofs
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getElementValues(std::string component, element_no_t elementNo, std::array<double,BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::nDofsPerElement()> &values)
{
  this->component_[component].getElementValues(elementNo, values);
}
*/
/*
//! get the values corresponding to all element-local dofs for all components
template<int D, typename BasisFunctionType>
template<int nComponents>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::nDofsPerElement()> &values)
{
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  std::vector<int> &dofGlobalNo = this->elementToDofMapping_->getElementDofs(elementNo);
  std::array<double,nComponents> resultVector;
  
  // transform global dof no.s to vector indices of first component
  for (int valueIndex = 0; valueIndex < nDofsPerElement; valueIndex++)
  {
    int valuesVectorIndex = dofGlobalNo[valueIndex]*nComponents;
    
    // create indices vector with values {0,1,2,...,nComponents-1}
    std::array<int,nComponents> indices;
    for(int i=0; i<nComponents; i++)
      indices[i] = valuesVectorIndex + i;
    
    // get values and assign them to result values vector
    VecGetValues(*this->values_, nComponents, indices.data(), resultVector.data());
    values[valueIndex] = resultVector;
  }
}
*/
//! for a specific component, get a single value from global dof no.
template<int D, typename BasisFunctionType>
double FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getValue(std::string component, node_no_t dofGlobalNo)
{
  assert(this->component_.find(component) != this->component_.end());
  return this->component_[component].getValue(dofGlobalNo); 
}

//! get a single value from global dof no. for all components
template<int D, typename BasisFunctionType>
template<int nComponents>
std::array<double,nComponents> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getValue(node_no_t dofGlobalNo)
{
  std::array<double,nComponents> resultVector;
  
  // transform global dof no.s to vector indices of first component
  int valuesVectorIndex = dofGlobalNo*nComponents;
  
  // create indices vector with values {0,1,2,...,nComponents-1}
  std::array<int,nComponents> indices;
  for(int i=0; i<nComponents; i++)
    indices[i] = valuesVectorIndex + i;
  
  // get values and assign them to result values vector
  VecGetValues(*this->values_, nComponents, indices.data(), resultVector.data());
  return resultVector;
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
initializeValuesVector()
{
  // initialize the PETSc vector that contains all the value entries
  // get number of entries
  this->nEntries_ = 0;
  // loop over components
  for (auto &component : this->component_)
  {
    this->nEntries_ += component.second.nDofs();
    VLOG(1) << "  component " << component.first << " has " << component.second.nDofs() << " dofs";
  }
  VLOG(1) << "total entries: " << this->nEntries_;
  
  // create vector
  this->values_ = std::make_shared<Vec>();
  
  PetscErrorCode ierr;
  // initialize PETSc vector object
  ierr = VecCreate(PETSC_COMM_WORLD, &*this->values_);  CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) *this->values_, this->name_.c_str()); CHKERRV(ierr);
  
  // initialize size of vector
  ierr = VecSetSizes(*this->values_, PETSC_DECIDE, this->nEntries_); CHKERRV(ierr);
  
  // set sparsity type and other options
  ierr = VecSetFromOptions(*this->values_);  CHKERRV(ierr);
      
  // set vector for all components
  for (auto &component : this->component_)
  {
    component.second.setValuesVector(this->values_);
  }
}

template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nComponents() const
{
  return this->component_.size(); 
}

template<int D, typename BasisFunctionType>
std::vector<std::string> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
componentNames() const
{
  std::vector<std::string> result;
  result.reserve(this->nComponents());
  for (auto componentIndex : this->component_)
  {   
    result.push_back(componentIndex.first);
  }
  return result;
}
  
template<int D, typename BasisFunctionType>
std::size_t FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nEntries() const
{
  return this->nEntries_; 
}

template<int D, typename BasisFunctionType>
dof_no_t FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nDofs() const
{
  return this->nEntries_ / this->nComponents(); 
}
  
template<int D, typename BasisFunctionType>
node_no_t FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nNodes() const
{
  return nodeToDofMapping_->nNodes();
}
  
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
initializeComponents(std::vector<std::string> &componentNames, std::string exfileBasisRepresentation)
{
  // insert components
  int nComponents = componentNames.size();
  int componentIndex = 0;
  
  // create a new values vector for the new field variable
  for(auto &componentName : componentNames)
  {
    Component<BasisOnMeshType> component;
    component.initialize(this->values_, nComponents, componentIndex++, this->nElements_);   // note: this->values_ may be nullptr but is updated by initializeValuesVector
    component.setName(componentName, exfileBasisRepresentation);
    component.setDofMappings(this->elementToDofMapping_, this->nodeToDofMapping_);
    component.setExfileRepresentation(this->exfileRepresentation_);
    this->component_.insert(std::pair<std::string, Component<BasisOnMeshType>>(componentName, component));
  }
  
  // set values_ and nEntries_
  initializeValuesVector();
}

template<int D, typename BasisFunctionType>
template<typename FieldVariableType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
initializeFromFieldVariable(FieldVariableType &fieldVariable, std::string name, std::vector<std::string> componentNames)
{
  this->name_ = name;
  this->exfileNo_ = 0;
  this->nEntries_ = 0;
  this->isGeometryField_ = false;
  this->nElements_ = fieldVariable.nElements();
  this->exfileRepresentation_ = fieldVariable.exfileRepresentation();
  this->elementToDofMapping_ = fieldVariable.elementToDofMapping();
  this->elementToNodeMapping_ = fieldVariable.elementToNodeMapping();
  this->nodeToDofMapping_ = fieldVariable.nodeToDofMapping();
  this->mesh_ = fieldVariable.mesh();
  
  std::string exfileBasisRepresentation = fieldVariable.component().begin()->second.exfileBasisFunctionSpecification();
  
  initializeComponents(componentNames, exfileBasisRepresentation);
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
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
  this->nElements_ = elementToDofMapping->nElements();
  
  // remove duplicate exfile representations
  exfileRepresentation->unifyExfileElementRepresentations();
  
  this->exfileRepresentation_ = exfileRepresentation;
  this->elementToDofMapping_ = elementToDofMapping;
  this->elementToNodeMapping_ = elementToNodeMapping;
  this->nodeToDofMapping_ = nodeToDofMapping;
  
  // this->mesh still needs to be set by setMesh
  
  std::string exfileBasisRepresentation = BasisFunction::getBasisRepresentationString<D,BasisFunctionType>();
  
  initializeComponents(componentNames, exfileBasisRepresentation);
  
  
  LOG(DEBUG) << "FieldVariable nDofs: " << this->nDofs();
}
  
template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getNumberScaleFactors(element_no_t globalElementNo) const
{
  
  //! return the node numbers and scale factors of the element
  return elementToNodeMapping_->getElement(globalElementNo).scaleFactors.size();
  
  //return component_.begin()->second.getNumberScaleFactors(globalNodeNo);
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
outputHeaderExelem(std::ostream &stream, element_no_t currentElementGlobalNo, int fieldVariableNo)
{
  // set no of field variable in ex file if it was specified
  if (fieldVariableNo != -1)
    exfileNo_ = fieldVariableNo+1;
  
  // output first line of header
  stream << " " << exfileNo_ << ") " << this->name_ << ", " << (isGeometryField_? "coordinate" : "field") 
    << ", rectangular cartesian, #Components=" << component_.size() << std::endl;
    
  // output headers of components
  for (auto &component : component_)
  {
    component.second.outputHeaderExelem(stream, currentElementGlobalNo);
  }
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
outputHeaderExnode(std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo)
{
  // set no of field variable in ex file if it was specified
  if (fieldVariableNo != -1)
    exfileNo_ = fieldVariableNo+1;
  
  // output first line of header
  stream << " " << exfileNo_ << ") " << this->name_ << ", " << (isGeometryField_? "coordinate" : "field") 
    << ", rectangular cartesian, #Components=" << component_.size() << std::endl;
    
  // output headers of components
  for (auto &component : component_)
  {
    component.second.outputHeaderExnode(stream, currentNodeGlobalNo, valueIndex);
  }
}

template<int D, typename BasisFunctionType>
bool FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
haveSameExfileRepresentation(element_no_t element1, element_no_t element2)
{
  // loop over components
  for (auto &component : component_)
  {
    if(!component.second.exfileRepresentation()->haveSameExfileRepresentation(element1, element2))
    {
      LOG(DEBUG) << "  component " << component.first << " has different exfileRepr for elements " << element1 << " and " << element2;
      return false;
    }
  }
  return true;
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
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
          double value = component.second.getValue(nodeDofs[dofIdx]);
          
          double scaleFactor2 = scaleFactors[dofIdx];
          assert(scaleFactor - scaleFactor2 < 1e-16);
          
          LOG(DEBUG) << "scaleFactor: " << scaleFactor << "," << scaleFactor2;
          
          // multiply value with scale factor
          value *= scaleFactor;
          
          nodeValues[dofIdx] = value;
        }
        
        // set updated values
        component.second.setValuesForNode(nodeGlobalNo, nodeValues.begin());
      }
    }
  }
}

template<int D, typename BasisFunctionType>
bool FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
isGeometryField()
{
  return isGeometryField_; 
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
output(std::ostream &stream) const
{
  stream << "\"" << this->name_ << "\", nEntries: " << nEntries_ << ", nElements: " << nElements_ 
    << ", isGeometryField: " << std::boolalpha << isGeometryField_ << std::endl
    << "  components:" << std::endl;
  for (auto &component : component_)
  {
    stream << component.second;
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

template<int D, typename BasisFunctionType>
std::ostream &operator<<(std::ostream &stream, const FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> &rhs)
{
  rhs.output(stream);
  return stream;
}
};  // namespace

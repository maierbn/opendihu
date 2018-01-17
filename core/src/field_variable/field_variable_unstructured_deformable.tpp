#include "field_variable/field_variable_unstructured_deformable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <map>
#include <fstream>

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
  bool newComponentStarts = false;
  int nComponents = 0;
  int componentNo = 0;
  int nLinesOfNodesFollow = 0;
  
  std::string componentName = "";
  std::string componentContent = "";
  
  // loop over content line-wise
  int pos = 0;
  while(pos < content.size())
  {
    // extract next line
    int posNewline = content.find("\n",pos);
    std::string line = content.substr(pos, posNewline);
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
      
      name_ = extractUntil(line, ",");
      std::string type = extractUntil(line, ",");
      isGeometryField_ = type.find("coordinate") != std::string::npos;
      nComponents = getNumberAfterString(line, "#Component<BasisOnMeshType>s=");
      newComponentStarts = true;
      componentNo = 0;
      continue;
    }
    
    // new block of component definition starts
    if (newComponentStarts)
    {
      // finish previous component
      if (componentNo != 0)
      {
        component_[componentName].initialize(values_, nComponents, componentNo, nElements_);
        component_[componentName].parseFromExelemFile(componentContent);
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
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
parseElementFromExelemFile(std::string content)
{
  for (typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter = component_.begin(); iter != component_.end(); iter++)
  {
    Component<BasisOnMeshType> &component = (*iter)->second;
    component.parseElementFromExelemFile(content);
  }
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
setNumberElements(element_idx_t nElements)
{
  this->nElements_ = nElements;
}

template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nElements() const
{
  return this->nElements_;
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
unifyMappings(ElementToNodeMapping &elementToNodeMapping, const int nDofsPerNode)
{
  // unify exfileRepresentation object
  // check if exfileRepresentation is equal for all components
  typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter = component_.begin();
  if (iter != component_.end())
  { 
    Component<BasisOnMeshType> &firstComponent = (*iter)->second;
    std::string firstComponentKey = (*iter)->first;
    
    this->exfileRepresentation_ = firstComponent.exfileRepresentation();
    
    // loop over further component after the first component
    iter++;
    for (; iter != component_.end(); iter++)
    {
      Component<BasisOnMeshType> &component = (*iter)->second;
      
      if (*this->exfileRepresentation_ == *component.exfileRepresentation())
      {
        component.setExfileRepresentation(this->exfileRepresentation_);
      }
      else 
      {
        LOG(ERROR) << "The components \"" << firstComponentKey << "\" and \"" << (*iter)->first 
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
  // loop over own components
  for (typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter = component_.begin(); iter != component_.end(); iter++)
  {
    Component<BasisOnMeshType> &component = (*iter)->second;
    
    // loop over components of the other fieldVariable
    for (typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter2 = fieldVariable.component_.begin(); iter2 != fieldVariable.component_.end(); iter2++)
    {
      Component<BasisOnMeshType> &component2 = (*iter2)->second;
      
      if (component.exfileRepresentation() == component2.exfileRepresentation())
      {
        component2.setExfileRepresentation(component.exfileRepresentation());
      }
      if (component.elementToDofMapping() == component2.elementToDofMapping())
      {
        component2.setElementToDofMapping(component.elementToDofMapping());
      }
    }
  }
  
}
  
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
createElementToDofMapping(ElementToNodeMapping &elementToNodeMapping, const int nDofsPerNode)
{
  
  // create the element to dof mapping from exfileRepresentation and element to node mapping
  this->nodeToDofMapping_ = this->elementToDofMapping_->setup(this->exfileRepresentation_, elementToNodeMapping, nDofsPerNode);
  
  // set the element to dof and node to dof mapping at each component
  for (typename std::map<std::string, Component<BasisOnMeshType>>::iterator iter = component_.begin(); iter != component_.end(); iter++)
  {
    Component<BasisOnMeshType> &component = (*iter)->second;
    component.setDofMappings(this->elementToDofMapping_, this->nodeToDofMapping_);
  }
}
  
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
parseFromExnodeFile(std::string content)
{
  int nFields = 0;
  int fieldNo = 0;
  std::stringstream nextFieldNo;
  std::map<std::string, int> componentDataBlockStartIndex;
  std::vector<double> values;
  int nodeNo = 0;
  
  enum {
    nothing,
    fieldsFollow,
    componentsFollow,
    valuesFollow,
  } lineType = nothing;
  
  // loop over file content line-wise
  int pos = 0;
  while(pos < content.size())
  {
    // extract next line
    int posNewline = content.find("\n",pos);
    std::string line = content.substr(pos, posNewline);
    if (posNewline == std::string::npos)
      pos = content.size();
    else
      pos = posNewline+1;
    
    bool finishPreviousNodeDataBlock = false;
    
    // line with number of fields, e.g. "#Fields=2", this indicates starting field description block
    if (line.find("#Fields=") != std::string::npos)
    {
      // finish previous values block
      if (lineType == valuesFollow)
      {
        // finish previous node
        finishPreviousNodeDataBlock = true;
      }
      
      // parse number of fields in current field description block
      nFields = getNumberAfterString(line, "#Fields=");
      
      // prepare string to look for for next field description
      fieldNo = 0;
      nextFieldNo.str("");
      nextFieldNo << fieldNo+1 << ")";
      lineType = fieldsFollow;
    }
    
    // line with field, e.g. "1) coordinates, coordinate, rectangular cartesian, #Component<BasisOnMeshType>s=3"
    if (line.find(nextFieldNo.str()) != std::string::npos)
    {
      // extract name of field
      extractUntil(line, ")");
      std::string name = extractUntil(line, ",");
      trim(name);
      
      // consider it if it matches the own name
      if (name == name_)
      {
        lineType = componentsFollow;
      }
      
      // prepare string to search for for the next field
      fieldNo++;
      nextFieldNo.str("");
      nextFieldNo << fieldNo+1 << ")";
      continue;
    }
    
    // line with component, e.g. "x.  Value index=1, #Derivatives=7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions=2"
    if (lineType == componentsFollow)
    {
      std::string componentName = extractUntil(line, ".");
      trim(componentName);
      
      if (component_.find(componentName) == component_.end())
      {
        LOG(WARNING) << "Component \"" << componentName << "\" of field variable \"" << name_ << "\" in exnode file is not present in exelem file.";
      }
      else
      {
        componentDataBlockStartIndex[componentName] = getNumberAfterString(line, "index=");
      }
    }
    
    // line with beginning of node, e.g. "Node: 73"
    if (line.find("Node:") != std::string::npos)
    {
      nodeNo = getNumberAfterString(line, "Node:");
      if (lineType == valuesFollow)
      {
        // finish previous node
        finishPreviousNodeDataBlock = true;
      }
      
      lineType = valuesFollow;
    }
      
    // if a node values block ended prior to the current line, store the collected values for the component
    if (finishPreviousNodeDataBlock)
    {
      // loop over components of current field variable
      for (auto &pair : componentDataBlockStartIndex)
      {
        // get data block index of component
        std::string componentName = pair.first;
        int index = pair.second;
        
        // store data at component
        component_[componentName].setValues(nodeNo, values.begin()+index);
      }
      // clear temporary data vector
      values.clear();
    }
    
    // values 
    if (lineType == valuesFollow)
    {
      trim(line);
      while(!line.empty())
      {
        values.push_back(atof(line.c_str()));
        extractUntil(line, " ");
        trim(line);
      }
    }
  }   // while
  
  
  if (lineType == valuesFollow)
  {
    // loop over components of current field variable
    for (auto &pair : componentDataBlockStartIndex)
    {
      // get data block index of component
      std::string componentName = pair.first;
      int index = pair.second;
      
      // store data at component
      component_[componentName].setValues(nodeNo, values.begin()+index);
    }
    // clear temporary data vector
    values.clear();
  }
}

template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getDofNo(element_idx_t elementNo, int dofIndex) const
{
  return this->component_.begin()->second.getDofNo(elementNo, dofIndex);
}

template<int D, typename BasisFunctionType>
Component<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> &FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
component(std::string key) const
{
  return this->component_[key]; 
}
  
template<int D, typename BasisFunctionType>
std::map<std::string, Component<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>> &FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
component() const
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
Vec &FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
values()
{
  return *this->values_; 
}

//! for a specific component, get values from their global dof no.s
template<int D, typename BasisFunctionType>
template<int N>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getValues(std::string component, std::array<int,N> dofGlobalNo, std::array<double,N> &values)
{
  this->component_[component].getValues(dofGlobalNo, values);
}
/*
template<int D, typename BasisFunctionType>
template<int N, int nComponents>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getValues(std::array<int,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
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
getElementValues(std::string component, element_idx_t elementNo, std::array<double,BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::nDofsPerElement()> &values)
{
  this->component_[component].getElementValues(elementNo, values);
}
*/
/*
//! get the values corresponding to all element-local dofs for all components
template<int D, typename BasisFunctionType>
template<int nComponents>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getElementValues(element_idx_t elementNo, std::array<std::array<double,nComponents>,BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::nDofsPerElement()> &values)
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
getValue(std::string component, node_idx_t dofGlobalNo)
{
  return this->component_[component].getValue(dofGlobalNo); 
}

//! get a single value from global dof no. for all components
template<int D, typename BasisFunctionType>
template<int nComponents>
std::array<double,nComponents> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getValue(node_idx_t dofGlobalNo)
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
    this->nEntries_ += component.second->nDofs();
  }
  
  // create vector
  this->values_ = std::make_shared<Vec>();
  
  PetscErrorCode ierr;
  // initialize PETSc vector object
  ierr = VecCreate(PETSC_COMM_WORLD, &*this->values_);  CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) *this->values_, this->name_); CHKERRV(ierr);
  
  // initialize size of vector
  ierr = VecSetSizes(*this->values_, PETSC_DECIDE, this->nEntries_); CHKERRV(ierr);
  
  // set sparsity type and other options
  ierr = VecSetFromOptions(*this->values_);  CHKERRV(ierr);
      
}

template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nComponents() const
{
  return this->component_.size(); 
}
  
template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nEntries() const
{
  return this->nEntries_; 
}

template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nDofs() const
{
  return this->nEntries_ / this->nComponents(); 
}
  
template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
nNodes() const
{
  return nodeToDofMapping_->nNodes();
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
  this->mesh_ = fieldVariable.mesh();
  
  // insert components
  int nComponents = componentNames.size();
  int componentIndex = 0;
  
  std::string exfileBasisRepresentation = fieldVariable->components().begin()->exfileBasisRepresentation();
  
  // create a new values vector for the new field variable
  PetscErrorCode ierr;
  
  for(auto &componentName : componentNames)
  {
    Component<BasisOnMeshType> component;
    component.initialize(this->values_, nComponents, componentIndex++, this->nElements_);
    component.setName(componentName, exfileBasisRepresentation);
    component.setElementToDofMapping(fieldVariable.elementToDofMapping());
    component.setExfileRepresentation(fieldVariable.exfileRepresentation());
    this->component_.insert(std::pair<std::string, Component<BasisOnMeshType>>(componentName, component));
  }
  
  // set values_ and nEntries_
  initializeValuesVector();
}
  
template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
getNumberScaleFactors() const
{
  return component_.begin()->second.getNumberScaleFactors();
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
outputHeaderExelem(std::ostream &stream, element_idx_t currentElementGlobalNo)
{
  stream << " " << exfileNo_ << ") " << name_ << ", " << (isGeometryField_? "coordinate" : "field") 
    << ", rectangular cartesian, #Components=" << component_.size() << std::endl;
    
  for (auto &component : component_)
  {
    component.second.outputHeaderExelem(stream, currentElementGlobalNo);
  }
}

template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
outputHeaderExnode(std::ostream &stream, node_idx_t currentNodeGlobalNo, int &valueIndex)
{
  stream << " " << exfileNo_ << ") " << name_ << ", " << (isGeometryField_? "coordinate" : "field") 
    << ", rectangular cartesian, #Components=" << component_.size() << std::endl;
    
  for (auto &component : component_)
  {
    component.second.outputHeaderExnode(stream, currentNodeGlobalNo, valueIndex);
  }
}

template<int D, typename BasisFunctionType>
bool FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
haveSameExfileRepresentation(element_idx_t element1, element_idx_t element2)
{
  // loop over components
  for (auto &component : component_)
  {
    if(!component.second.exfileRepresentation().haveSameExfileRepresentation(element1, element2))
    {
      return false;
    }
  }
  return true;
}

};  // namespace
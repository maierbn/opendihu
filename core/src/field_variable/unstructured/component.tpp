#include "field_variable/unstructured/component.h"

#include <petscvec.h>
#include "utility/string_utility.h"

namespace FieldVariable
{

using namespace StringUtility;
  
template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
initialize(std::shared_ptr<Vec> values, int nComponents, int componentIndex, int nElements)
{
  values_ = values;
  nComponents_ = nComponents;
  componentIndex_ = componentIndex;
  nElements_ = nElements;
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
setValuesVector(std::shared_ptr<Vec> values)
{
  assert(values);
  values_ = values;
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
parseHeaderFromExelemFile(std::string content)
{
  VLOG(1) << "parseHeaderFromExelemFile";
  int firstLine = true;
  std::string exfileRepresentationContent;
  
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
    
    if (firstLine) 
    {
      name_ = extractUntil(line, ".");
      trim(name_);
      VLOG(1) << "component " << name_;
      exfileBasisFunctionSpecification_ = extractUntil(line, ",");
      trim(exfileBasisFunctionSpecification_);
      firstLine = false;
    }
    else 
    {
      exfileRepresentationContent += line+"\n";
    }
  }
  
  if (!exfileRepresentation_)
  {
    exfileRepresentation_ = std::make_shared<ExfileRepresentation>();
    exfileRepresentation_->setNumberElements(nElements_);
  }
  exfileRepresentation_->parseHeaderFromExelemFile(exfileRepresentationContent);
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
parseElementFromExelemFile(std::string content)
{
  VLOG(1) << "component " << name_ << ", parseElementFromExelemFile";
  exfileRepresentation_->parseElementFromExelemFile(content);
  //VLOG(1) << "   after parseElementFromExelemFile: " << *exfileRepresentation_;
}

template<typename BasisOnMeshType>
std::shared_ptr<ExfileRepresentation> Component<BasisOnMeshType>::
exfileRepresentation()
{
  return exfileRepresentation_;
}

template<typename BasisOnMeshType>
std::shared_ptr<ElementToDofMapping> Component<BasisOnMeshType>::
elementToDofMapping()
{
  return elementToDofMapping_;
}

template<typename BasisOnMeshType>
std::shared_ptr<NodeToDofMapping> Component<BasisOnMeshType>::
nodeToDofMapping()
{
  return nodeToDofMapping_;
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
setExfileRepresentation(std::shared_ptr<ExfileRepresentation> exfileRepresentation)
{
  exfileRepresentation_ = exfileRepresentation;
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
setName(std::string name, std::string exfileBasisFunctionSpecification)
{
  this->name_ = name;
  this->exfileBasisFunctionSpecification_ = exfileBasisFunctionSpecification;
}

template<typename BasisOnMeshType>
std::string Component<BasisOnMeshType>::
exfileBasisFunctionSpecification() const
{
  return this->exfileBasisFunctionSpecification_;
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
setDofMappings(std::shared_ptr<ElementToDofMapping> elementToDofMapping, std::shared_ptr<NodeToDofMapping> nodeToDofMapping)
{
  elementToDofMapping_ = elementToDofMapping;
  nodeToDofMapping_ = nodeToDofMapping;
}
  
template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
setNodeValues(node_no_t nodeGlobalNo, std::vector<double>::iterator valuesBegin)
{
  NodeToDofMapping::NodeDofInformation &nodeDofInformation = nodeToDofMapping_->getNodeDofInformation(nodeGlobalNo);
  
  // loop over dofs
  for(int dofIndex = 0; dofIndex < nodeDofInformation.dofs.size(); dofIndex++, valuesBegin++)
  {
    dof_no_t dofGlobalNo = nodeDofInformation.dofs[dofIndex];
    double value = *valuesBegin;
    std::size_t vectorIndex = dofGlobalNo*nComponents_ + componentIndex_;
    
    VLOG(2) << " component " << name_ << ", set value: " << value << " at nodeGlobalNo: " << nodeGlobalNo 
      << ", componentIndex: " << componentIndex_ << ", nComponents: " << nComponents_ 
      << ", dofGlobalNo: " <<dofGlobalNo << " -> set at vectorIndex: " << vectorIndex;
    
    VecSetValue(*values_, vectorIndex, value, INSERT_VALUES);
  }
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
setNodeValuesFromBlock(node_no_t nodeGlobalNo, std::vector<double>::iterator valuesBegin)
{
  NodeToDofMapping::NodeDofInformation &nodeDofInformation = nodeToDofMapping_->getNodeDofInformation(nodeGlobalNo);
  
  std::vector<int> indices;
  std::vector<double> values;
  
  const unsigned int nVersions = nodeDofInformation.elementsOfVersion.size();
  indices.reserve(BasisOnMeshType::nDofsPerNode()*nVersions);
  values.reserve(BasisOnMeshType::nDofsPerNode()*nVersions);
  
  VLOG(2) << "setNodeValuesFromBlock(" << nodeGlobalNo << "), block: ";
  for (std::vector<double>::iterator iter = valuesBegin; iter < valuesBegin + BasisOnMeshType::nDofsPerNode()*nVersions; iter++)
    VLOG(2) << "  " << *iter;
  
  VLOG(2) << " reserve " << BasisOnMeshType::nDofsPerNode() << "*" << nVersions << " entries";
  // loop over versions 
  for (unsigned int versionIdx = 0; versionIdx < nVersions; versionIdx++)
  {
    VLOG(2) << "version " << versionIdx;
   
    // loop over elements of that version 
    for (auto element : nodeDofInformation.elementsOfVersion[versionIdx])
    {
      element_no_t elementGlobalNo = element.elementGlobalNo;
      node_no_t nodeIdx = element.nodeIdx;
      std::shared_ptr<ExfileElementRepresentation> exfileElement = exfileRepresentation_->getExfileElementRepresentation(elementGlobalNo);
     
      VLOG(2) << " element " << elementGlobalNo << " nodeIdx " << nodeIdx << ", node has the following valueIndices: " << exfileElement->getNode(nodeIdx).valueIndices;
      
      // loop over block indices
      for (unsigned int i = 0; i < exfileElement->getNode(nodeIdx).valueIndices.size(); i++)
      {
        int blockIndex = exfileElement->getNode(nodeIdx).valueIndices[i];
        VLOG(2) << "  block index " << blockIndex;
        
        double value = *(valuesBegin + blockIndex);
        
        VLOG(2) << "  value " << value;
        values.push_back(value);
        
        dof_no_t dofGlobalNo = nodeDofInformation.dofs[blockIndex];
        std::size_t vectorIndex = dofGlobalNo*nComponents_ + componentIndex_;
        VLOG(2) << "  dofGlobalNo = " << dofGlobalNo << ", vectorIndex = " << vectorIndex;
        
        indices.push_back(vectorIndex);
      }
    }
  }
  
  assert (values_);
  VLOG(2) << " set " << indices.size() << " values at indices " << indices;
  VLOG(2) << "                                " << values;
  VecSetValues(*values_, indices.size(), indices.data(), values.data(), INSERT_VALUES);
}

template<typename BasisOnMeshType>
dof_no_t Component<BasisOnMeshType>::nDofs() const
{
  return this->elementToDofMapping_->nDofs();
}

template<typename BasisOnMeshType>
element_no_t Component<BasisOnMeshType>::nElements() const
{
  return this->nElements_;
}

template<typename BasisOnMeshType>
dof_no_t Component<BasisOnMeshType>::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  return elementToDofMapping_->getElementDofs(elementNo)[dofIndex];
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
getValues(std::vector<double> &values, bool onlyNodalValues)
{
  const dof_no_t nDofs = this->nDofs();
 
  // set stride to 2 if Hermite, else to 1  
  const int stride = (onlyNodalValues && std::is_same<typename BasisOnMeshType::BasisFunction, BasisFunction::Hermite>::value ? 2 : 1);
  
  // create vector indices for all dofs
  dof_no_t nValues = nDofs;
  
  // if Hermite and only values at nodes should be retrieved, the number of values is half the number of dofs
  if (onlyNodalValues)
    // if the basis function is Hermite
    if (std::is_same<typename BasisOnMeshType::BasisFunction, BasisFunction::Hermite>::value)
      nValues = nDofs / 2;
    
  std::vector<int> indices(nValues,0);
  dof_no_t indexNo = 0;
  for (dof_no_t dofGlobalNo=0; dofGlobalNo<nDofs; dofGlobalNo+=stride)
  {
    assert(indexNo < nValues);
    indices[indexNo++] = dofGlobalNo*this->nComponents_ + componentIndex_;
  }
  
  values.resize(nValues);
  assert (values_);
  VecGetValues(*values_, nValues, indices.data(), values.data());
}

template<typename BasisOnMeshType>
template<int N>
void Component<BasisOnMeshType>::
getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values)
{
  // transform global dof no.s to vector indices
  for (auto &index : dofGlobalNo)
  {
    index = index*nComponents_ + componentIndex_;
  }
 
  assert (values_);
  VecGetValues(*values_, N, dofGlobalNo.data(), values.data());
}
  
template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
getValues(std::vector<dof_no_t> dofGlobalNo, std::vector<double> &values)
{
  const int nValues = dofGlobalNo.size();
  
  // transform global dof no.s to vector indices
  for (auto &index : dofGlobalNo)
  {
    index = index*nComponents_ + componentIndex_;
  }
 
 VLOG(1) << "Component getValues, " << nValues << " values, componentIndex=" << componentIndex_ 
   << ", nComponents= " << nComponents_ << " indices: " << dofGlobalNo;
   
  assert (values_);
  
  std::size_t previousSize = values.size();
  VLOG(1) << "previousSize: " << previousSize;
  values.resize(previousSize+nValues);
  VLOG(1) << "new size: " << values.size();
  VecGetValues(*values_, nValues, (PetscInt*)dofGlobalNo.data(), values.data()+previousSize);
  
  
  VLOG(1) << "retrieved values: " << values;
}
  
template<typename BasisOnMeshType>
double Component<BasisOnMeshType>::
getValue(node_no_t dofGlobalNo)
{
  double value;
  std::array<int,1> indices{(int)(dofGlobalNo*nComponents_ + componentIndex_)};
  
  assert (values_);
  VecGetValues(*values_, 1, indices.data(), &value);
  
  VLOG(2) << "component " << this->name_<<", getValue for dof " << dofGlobalNo 
    << " componentIndex: " << componentIndex_ << ", nComponents: " << nComponents_ << ", vectorIndex: " << indices[0] 
    << ", value: " << value;
  return value;
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
getElementValues(element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  const std::vector<dof_no_t> &elementDofs = elementToDofMapping_->getElementDofs(elementNo);
  
  std::array<PetscInt, BasisOnMeshType::nDofsPerElement()> indices;
  
  // transform global dof no.s to vector indices
  for (int dofIndex = 0; dofIndex < BasisOnMeshType::nDofsPerElement(); dofIndex++)
  {
    indices[dofIndex] = elementDofs[dofIndex] * nComponents_ + componentIndex_;
  }
  
  assert (values_);
  VecGetValues(*values_, BasisOnMeshType::nDofsPerElement(), indices.data(), values.data());
}

template<typename BasisOnMeshType>
int Component<BasisOnMeshType>::
getNumberScaleFactors(node_no_t nodeGlobalNo) const
{
  return nodeToDofMapping_->getNodeScaleFactors(nodeGlobalNo).size();
}
  
template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo)
{
  file << " " << name_ << ". " << exfileBasisFunctionSpecification_ << ", no modify, standard node based." << std::endl
    << "   #Nodes=" << BasisOnMeshType::nNodesPerElement() << std::endl;
  exfileRepresentation_->getExfileElementRepresentation(currentElementGlobalNo)->outputHeaderExelem(file);
}
  
template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex)
{
  int nDerivativesPerNode = BasisOnMeshType::nDofsPerNode() - 1;
  std::string derivativeSpecifiers = "";
  if (nDerivativesPerNode == 1)
    derivativeSpecifiers = " (d/ds1)";
  else if (nDerivativesPerNode == 3)
    derivativeSpecifiers = " (d/ds1,d/ds2,d2/ds1ds2)";
  else if (nDerivativesPerNode == 7)
    derivativeSpecifiers = " (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3)";

  int nVersions = nodeToDofMapping_->nVersions(currentNodeGlobalNo);
  
  file << "  " << name_ << ".  Value index=" << valueIndex+1 << ", #Derivatives=" << nDerivativesPerNode << derivativeSpecifiers 
    << ", #Versions=" << nVersions << std::endl;
  
  valueIndex += BasisOnMeshType::nDofsPerNode() * nVersions;
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
output(std::ostream &stream) const
{
  stream << "\"" << name_ << "\", componentIndex_: " << componentIndex_ << "/" << nComponents_ << ", nElements: " << nElements_ 
    << ", exfileBasisFunctionSpecification_: " << exfileBasisFunctionSpecification_ << std::endl
    << ", exfileRepresentation_=" << exfileRepresentation_ << ", nodeToDofMapping_=" << nodeToDofMapping_ 
    << ", elementToDofMapping_=" << elementToDofMapping_ << std::endl << "exfileRepresentation: " << *exfileRepresentation_;
}

template<typename BasisOnMeshType>
std::ostream &operator<<(std::ostream &stream, const Component<BasisOnMeshType> &rhs)
{
  rhs.output(stream);
  return stream; 
}

};

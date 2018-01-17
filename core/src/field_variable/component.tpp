#include "field_variable/component.h"

#include <petscvec.h>
#include "utility/string_utility.h"

namespace FieldVariable
{

using namespace StringUtility;
  
template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
initialize(std::shared_ptr<Vec> &values, int nComponents, int componentIndex, int nElements)
{
  values_ = values;
  nComponents_ = nComponents;
  componentIndex_ = componentIndex;
  nElements_ = nElements;
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
parseHeaderFromExelemFile(std::string content)
{
  int firstLine = true;
  std::string exfileRepresentationContent;
  
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
    
    if (firstLine) 
    {
      name_ = extractUntil(line, ".");
      trim(name_);
      exfileBasisFunctionSpecification_ = extractUntil(line, ",");
      trim(exfileBasisFunctionSpecification_);
      firstLine = false;
    }
    else 
    {
      exfileRepresentationContent += line+"\n";
    }
  }
  exfileRepresentation_ = std::make_shared<ExfileRepresentation>();
  exfileRepresentation_->setNumberElements(nElements_);
  exfileRepresentation_->parseHeaderFromExelemFile(exfileRepresentationContent);
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
parseElementFromExelemFile(std::string content)
{
  exfileRepresentation_->parseElementFromExelemFile(content);
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
  return elementToDofMapping;
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
setExfileRepresentation(std::shared_ptr< ExfileRepresentation > exfileRepresentation)
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
setValuesForNode(node_idx_t nodeGlobalNo, std::vector<double>::iterator valuesBegin)
{
  ElementToDofMapping::nodeDofInformation &nodeDofInformation = elementToDofMapping_->getNodeDofInformation(nodeGlobalNo);
  
  for(int i = 0; i < nodeDofInformation.exfileValueIndices; i++)
  {
    int valueBlockIndex = nodeDofInformation.exfileValueIndices[i];
    int dofNo = nodeDofInformation.dofs[i];
    
    double value = *(valuesBegin+valueBlockIndex);
    
    int vectorIndex = dofNo*nComponents_ + componentIndex_;
    
    VecSetValue(*values_, vectorIndex, value, INSERT_VALUES);
  }
}

template<typename BasisOnMeshType>
int Component<BasisOnMeshType>::nDofs() const
{
  return this->elementToDofMapping_->nDofs();
}

template<typename BasisOnMeshType>
int Component<BasisOnMeshType>::nElements() const
{
  return this->nElements_;
}

template<typename BasisOnMeshType>
int Component<BasisOnMeshType>::
getDofNo(element_idx_t elementNo, int dofIndex) const
{
  return elementToDofMapping_->getElementDofs(elementNo)[dofIndex];
}

template<typename BasisOnMeshType,int N>
void Component<BasisOnMeshType>::
getValues(std::array<int,N> dofGlobalNo, std::array<double,N> &values)
{
  // transform global dof no.s to vector indices
  for (auto &index : dofGlobalNo)
  {
    index = index*nComponents_ + componentIndex_;
  }
 
  VecGetValues(values_, N, dofGlobalNo.data(), values.data());
}
  
template<typename BasisOnMeshType>
double Component<BasisOnMeshType>::
getValue(node_idx_t dofGlobalNo)
{
  double value;
  std::array<int,1> indices{dofGlobalNo*nComponents_ + componentIndex_};
  VecGetValues(values_, 1, indices.data(), &value);
}

template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
getElementValues(element_idx_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  std::vector<int> &dofs = elementToDofMapping_->getElementDofs(elementNo);
  
  // transform global dof no.s to vector indices
  for (auto &index : dofs)
  {
    index = index*nComponents_ + componentIndex_;
  }
  
  VecGetValues(values_, BasisOnMeshType::nDofsPerElement(), dofs.data(), values.data());
}

template<typename BasisOnMeshType>
int Component<BasisOnMeshType>::
getNumberScaleFactors(node_idx_t nodeGlobalNo)
{
  return elementToDofMapping_->getNodeScaleFactors(nodeGlobalNo).size();
}
  
template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
outputHeaderExelem(std::ofstream file, element_idx_t currentElementGlobalNo)
{
  file << " " << name_ << ". " << exfileBasisFunctionSpecification_ << ", no modify, standard node based." << std::endl
    << " #Nodes=" << BasisOnMeshType::nNodesPerElement() << std::endl;
  file << exfileRepresentation_->getExfileElementRepresentation(currentElementGlobalNo).outputHeaderExelemFile(file);
}
  
template<typename BasisOnMeshType>
void Component<BasisOnMeshType>::
outputHeaderExnode(std::ostream file, node_idx_t currentNodeGlobalNo, int &valueIndex)
{
  int nDerivativesPerNode = BasisOnMeshType::nDofsPerNode() - 1;
  std::string derivativeSpecifiers = "";
  if (nDerivativesPerNode == 1)
    derivativeSpecifiers = " (d/ds1)";
  else if (nDerivativesPerNode == 3)
    derivativeSpecifiers = " (d/ds1,d/ds2,d2/ds1ds2)";
  else if (nDerivativesPerNode == 7)
    derivativeSpecifiers = " (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3)";

  int nVersions = nodeToDofMapping_->getNumberVersions(currentNodeGlobalNo, BasisOnMeshType::nNodesPerElement());
  
  file << " " << name_ << ".  Value index=" << valueIndex << ", #Derivatives=" << nDerivativesPerNode << derivativeSpecifiers 
    << ", #Versions=" << nVersions << std::endl;
  
  valueIndex += BasisOnMeshType::nDofsPerNode() * nVersions;
}

};
#include "basis_on_mesh/03_basis_on_mesh_dofs.h"

#include "easylogging++.h"
#include "utility/string_utility.h"

#include "basis_function/basis_function.h"
#include "field_variable/factory.h"

#include <iostream>
#include <fstream>

namespace BasisOnMesh
{
 
template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
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
    this->nDofs_ = elementToDofMapping->nDofs();
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
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
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
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
outputExelemFile(std::ostream &file)
{ 
  file << " Group name: Region" << std::endl
    << " Shape. Dimension=" << D << ", " << StringUtility::multiply<D>("line") << std::endl;
   
  bool outputHeader = true;
   
  for(element_no_t currentElementGlobalNo = 0; currentElementGlobalNo < nElements(); currentElementGlobalNo++)
  {
    int nScaleFactors = fieldVariable_.begin()->second->getNumberScaleFactors(currentElementGlobalNo);
    
    if (nScaleFactors != 0)
    {
      file << " #Scale factor sets=1" << std::endl;
      file << " " << BasisFunction::getBasisRepresentationString<D,BasisFunctionType>();
      file << ", #Scale factors=" << nScaleFactors << std::endl;
    }
    
    file << " #Nodes=" << BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::nNodesPerElement() << std::endl
      << " #Fields=" << fieldVariable_.size() << std::endl;
      
    // check if a new header is necessary
    if (currentElementGlobalNo > 0)
    {
      outputHeader = false;
      for (auto &fieldVariable : fieldVariable_)
      {
        if (!fieldVariable->second->haveSameExfileRepresentation(currentElementGlobalNo-1, currentElementGlobalNo))
        {
          outputHeader = true;
          break;
        }
      }
    }
    
    // output header 
    if (outputHeader)
    {
      for (auto &fieldVariable : fieldVariable_)
      {
        fieldVariable->second->outputHeaderExelem(file, currentElementGlobalNo);
      }
    }
    
    elementToNodeMapping_->outputElementExelem(file, currentElementGlobalNo);
  }
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
outputExnodeFile(std::ostream &file)
{ 
  file << " Group name: Region" << std::endl;

  const int nNodesPerElement = BasisOnMeshDofs<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::nNodesPerElement();   
  bool outputHeader = true;
   
  for(node_no_t currentNodeGlobalNo = 0; currentNodeGlobalNo < nElements(); currentNodeGlobalNo++)
  {
    // check if a new header is necessary
    if (currentNodeGlobalNo > 0)
    {
      outputHeader = false;
      for (auto &fieldVariable : fieldVariable_)
      {
        int previousNumberVersions = fieldVariable->second->nodeToDofMapping()->getNumberVersions(currentNodeGlobalNo-1, nNodesPerElement);
        int currentNumberVersions = fieldVariable->second->nodeToDofMapping()->getNumberVersions(currentNodeGlobalNo, nNodesPerElement);
        
        if(previousNumberVersions != currentNumberVersions)
        {
          outputHeader = true;
          break;
        }
      }
    }
    
    // output header 
    if (outputHeader)
    {
      file << " #Fields=" << fieldVariable_.size() << std::endl;
      int valueIndex = 0;  // an index that runs over values of components of field variables and corresponds to the index in the values block for a node
      for (auto &fieldVariable : fieldVariable_)
      {
        fieldVariable->second->outputHeaderExnode(file, currentNodeGlobalNo, valueIndex);
      }
    }
    
    file << " Node: " << currentNodeGlobalNo << std::endl;
    
    // collect values of all field variables at the current node
    // get dofs
    std::vector<double> valuesAtNode;
    for (auto &fieldVariable : fieldVariable_)
    {
      std::vector<int> dofsAtNode;
      fieldVariable->second->mesh()->getNodeDofs(currentNodeGlobalNo, dofsAtNode);
      
      // loop over components
      for (int componentNo = 0; componentNo <fieldVariable->second.getNComponents(); componentNo++)
      {
        // loop over dofs 
        for (dof_no_t dofGlobalNo : dofsAtNode)
        {
          double value = fieldVariable->second->getValue(componentNo, dofGlobalNo);
          valuesAtNode.push_back(value);
        }
      }
    }
    
    StringUtility::outputValuesBlock(file, valuesAtNode.begin(), valuesAtNode.end(), 8);
  }
}

};  // namespace

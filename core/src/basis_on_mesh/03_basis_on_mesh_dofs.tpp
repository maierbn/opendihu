#include "basis_on_mesh/03_basis_on_mesh_dofs.h"

#include <cmath>
#include <array>
#include <string>

#include "easylogging++.h"
#include "utility/string_utility.h"

namespace BasisOnMesh
{
 
using namespace StringUtility;

// element-local dofIndex to global dofNo for 1D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getDofNo(element_idx_t elementNo, int dofIndex) const
{
  return BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>>::getDofNo(this->nElements_, elementNo, dofIndex);
}


template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getDofNo(std::array<element_idx_t, MeshType::dim()> nElements, element_idx_t elementNo, int dofIndex)
{
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0,1 2,3
  // averageNDofsPerElement: 
  // 1         2            2
  return BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNDofsPerElement() * elementNo + dofIndex;
}
  
 

// element-local dofIndex to global dofNo for 2D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getDofNo(element_idx_t elementNo, int dofIndex) const
{
  return BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>>::getDofNo(this->nElements_, elementNo, dofIndex);
}

// element-local dofIndex to global dofNo for 2D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getDofNo(std::array<element_idx_t, MeshType::dim()> nElements, element_idx_t elementNo, int dofIndex)
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

// element-local dofIndex to global dofNo for 3D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getDofNo(element_idx_t elementNo, int dofIndex) const
{
  return BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>>::getDofNo(this->nElements_, elementNo, dofIndex);
}
 
// element-local dofIndex to global dofNo for 3D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getDofNo(std::array<element_idx_t, MeshType::dim()> nElements, element_idx_t elementNo, int dofIndex)
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

// element-local nodeIndex to global nodeNo for 1D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeNo(element_idx_t elementNo, int nodeIndex) const
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
getNodeNo(element_idx_t elementNo, int nodeIndex) const
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      2 3
  // 0 1       0 1 2      0 1
  // nNodesPerElement:
  // 4         9          4
  
  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given
  const std::array<int,2> nElements{this->nElements(0), this->nElements(1)};
  
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
getNodeNo(element_idx_t elementNo, int nodeIndex) const
{
  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given
  const std::array<int,3> nElements{this->nElements(0), this->nElements(1), this->nElements(2)};
  
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
element_idx_t BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
nElements() const
{
  return this->nElements_;
}

template<int D,typename BasisFunctionType>
BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
BasisOnMeshDofs(PyObject *settings) :
  BasisOnMeshJacobian<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::BasisOnMeshJacobian(settings)
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
}


template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
initialize()
{
  LOG(DEBUG) << "   retrieve this pointer ";
  std::shared_ptr<BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> ptr = this->shared_from_this();
  
  assert(ptr != nullptr);
  
  LOG(DEBUG) << "   cast this pointer ";
  std::shared_ptr<BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> self = std::static_pointer_cast<BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>(ptr);
  
  assert(self != nullptr);
  this->geometry_->setMesh(self);
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
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
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
parseExelemFile(std::string exelemFilename)
{
  std::ifstream file_exelem(exelemFilename.c_str(), std::ios::in | std::ios::binary);
  if (!file_exelem.is_open())
  {
    LOG(WARNING) << "Could not open exelem file \"" << exelemFilename << "\" for reading.";
  }
   
  // find out number of elements in file
  this->nElements_ = 0;
  std::string line;
  while(!file_exelem.eof())
  {
    getline(file_exelem, line);
    if(file_exelem.eof())
      break;
    
    if (line.find("Element:") != std::string::npos)
    {
      int elementGlobalNo = getNumberAfterString(line, "Element:");
      this->nElements_ = std::max(this->nElements_, elementGlobalNo);
    }
  }
  this->elementToNodeMapping_.setNumberElements(this->nElements_);
  
  // read in dofs for each element
  int nFields = 0;
  int fieldNo = 0;
  int nComponents = 0;
  int nNodes = 0;
  int nodeNo = 0;
  int nValues = 0;
  
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
  
  while(!file_exelem.eof())
  {
    getline(file_exelem, line);
    if(file_exelem.eof())
      break;
    
    // check if line contains "Shape."
    if (line.find("Shape.") != std::string::npos && line.find("Dimension=") != std::string::npos)
    {
      int dimension = getNumberAfterString(line, "Dimension=");
      if (dimension == 3)
      {
        lineType = elementsFollow;
      }
      continue;
    }
    
    // check if line contains "#Fields="
    if(lineType == elementsFollow && line.find("#Fields=") != std::string::npos) 
    {
      // parse number of following fields
      nFields = getNumberAfterString(line, "#Fields=");
      
      // prepare string that starts the next field description block
      fieldNo = 0;
      nextFieldNo.str("");
      nextFieldNo << fieldNo+1 << ")";
      lineType = fieldsFollow;
      continue;
    }
    
    // if in field description block
    if (lineType == fieldsFollow)
    {
      // if the next field description block begins or a block with elements begins (i.e. the current field description block is finished)
      if(line.find(nextFieldNo.str()) != std::string::npos || line.find("Element:") != std::string::npos)
      {
        // finish previous fieldVariable component
        if (fieldNo != 0)
        {
          if (this->fieldVariable_.find(fieldName) == this->fieldVariable_.end())
          {
            this->fieldVariable_.insert(std::pair<std::string, 
                                        std::shared_ptr<FieldVariableType>>(fieldName, 
              std::make_shared<BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>()));
          }
          this->fieldVariable_[fieldName]->setNumberElements(this->nElements_);
          
          // parse whole field description block inside component object
          this->fieldVariable_[fieldName]->parseFromExelemFile(fieldContent);
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
        // if there were lines collected for the previous element block
        if (!elementContent.empty())
        {
          // parse whole element block into elmeentToNode mapping
          this->elementToNodeMapping_.parseElementFromExelemFile(elementContent);
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
    this->elementToNodeMapping_.parseElementFromExelemFile(elementContent);
    for(auto &fieldVariable : this->fieldVariable_)
    {
      fieldVariable.second->parseElementFromExelemFile(elementContent);
    }
  }
  
  // unify exfile representation variables in field variables such that there is only one shared for all components of a field variable
  for(auto &fieldVariable : this->fieldVariable_)
  {
    fieldVariable.second->unifyMappings(this->elementToNodeMapping_, this->nDofsPerNode());
  }
  
  // remove duplicate exfile representation and element to dof mapping entries that are shared among different field variables
  for(auto &fieldVariable : this->fieldVariable_)
  {
    for(auto &fieldVariable2 : this->fieldVariable_)
    {
      if (fieldVariable.first != fieldVariable2.first)
      {
        fieldVariable.second->unifyMappings(fieldVariable2.second);
      }
    }
  }
  
  // allocate common values vector for each field variable
  for(auto &fieldVariable : this->fieldVariable_)
  {
    fieldVariable.second->initializeValuesVector();
  }
}
  
template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
parseExnodeFile(std::string exnodeFilename)
{
  std::ifstream file_exnode(exnodeFilename.c_str(), std::ios::in | std::ios::binary);
  if (!file_exnode.is_open())
  {
    LOG(WARNING) << "Could not open exnode file \"" << exnodeFilename << "\" for reading.";
  }
   
  std::string content( (std::istreambuf_iterator<char>(file_exnode) ),
                       (std::istreambuf_iterator<char>()    ) );
  
  for(auto &fieldVariable : this->fieldVariable_)
  {
    fieldVariable.second->parseFromExnodeFile(content);
  }
}

template<int D,typename BasisFunctionType>
int BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
getDofNo(element_idx_t elementNo, int dofIndex) const
{
  if (this->fieldVariable_.find("geometry") == this->fieldVariable_.end())
    LOG(FATAL) << "Mesh contains no field variable \"geometry\". Use remap to create one!";
  
  return this->fieldVariable_["geometry"]->getDofNo(elementNo, dofIndex);
}

template<int D,typename BasisFunctionType>
int BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
getNodeNo(element_idx_t elementNo, int nodeIndex) const
{
  return this->elementToNodeMapping_->getElement(elementNo).globalNodeNo[nodeIndex];
}

};  // namespace
#include "basis_on_mesh/04_basis_on_mesh_data_unstructured.h"

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
BasisOnMeshDataUnstructured<D,BasisFunctionType>::
BasisOnMeshDataUnstructured(std::shared_ptr<Partition::Manager> partitionManager, PyObject *settings, bool noGeometryField) :
  BasisOnMeshPartition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::BasisOnMeshPartition(partitionManager),
  noGeometryField_(noGeometryField)
{
  LOG(TRACE) << "BasisOnMeshDataUnstructured constructor";

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
dof_no_t BasisOnMeshDataUnstructured<D,BasisFunctionType>::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  if (!this->geometryField_)
    LOG(FATAL) << "Mesh contains no field variable \"geometry\". Use remap to create one!";

  return this->geometryField_->getDofNo(elementNo, dofIndex);
}

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDataUnstructured<D,BasisFunctionType>::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  return this->elementToNodeMapping_->getElement(elementNo).nodeGlobalNo[nodeIndex];
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDataUnstructured<D,BasisFunctionType>::
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
dof_no_t BasisOnMeshDataUnstructured<D,BasisFunctionType>::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const 
{
  return this->geometryField_->nodeToDofMapping()->getNodeDofNo(nodeGlobalNo, dofIndex);
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDataUnstructured<D,BasisFunctionType>::
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
void BasisOnMeshDataUnstructured<D,BasisFunctionType>::
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
std::shared_ptr<FieldVariable::ElementToNodeMapping> BasisOnMeshDataUnstructured<D,BasisFunctionType>::
elementToNodeMapping()
{
  return elementToNodeMapping_;
}

template<int D,typename BasisFunctionType>
element_no_t BasisOnMeshDataUnstructured<D,BasisFunctionType>::
nElements() const
{
  return this->nElements_;
}
};  // namespace

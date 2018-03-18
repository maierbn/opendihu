#include "field_variable/00_field_variable_base.h"

namespace FieldVariable
{

template<typename BasisOnMeshType> 
FieldVariableBase<BasisOnMeshType>::
FieldVariableBase() : mesh_(nullptr)
{
}

template<typename BasisOnMeshType> 
std::shared_ptr<BasisOnMeshType> FieldVariableBase<BasisOnMeshType>::
mesh()
{
  return mesh_;
}

template<typename BasisOnMeshType> 
void FieldVariableBase<BasisOnMeshType>::
setMesh(std::shared_ptr<BasisOnMeshType> mesh)
{
  mesh_ = mesh;
}

template<typename BasisOnMeshType> 
std::string FieldVariableBase<BasisOnMeshType>::
name() const
{
  return this->name_;
}

} // namespace 

#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>

namespace FieldVariable
{
 
/** base class for a field variable that just stores the mesh the field variable is defined on
 */
template<typename BasisOnMeshType>
class FieldVariableBase
{
public:
  FieldVariableBase();
  
  //! return the mesh of this field variable
  std::shared_ptr<BasisOnMeshType> mesh();
  
  //! set the internal mesh
  void setMesh(std::shared_ptr<BasisOnMeshType> mesh);
  
  //! get the name of the field variable
  std::string name() const;
  
protected:
  std::shared_ptr<BasisOnMeshType> mesh_;  ///< the mesh for which the field variable is defined
  std::string name_;     ///< name of the field variable
  
};
 
}; // namespace 
#include "field_variable/field_variable_base.tpp"

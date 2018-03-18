#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>

#include "field_variable/unstructured/01_field_variable_data_unstructured_deformable.h"
#include "field_variable/unstructured/component.h"
#include "field_variable/unstructured/element_to_node_mapping.h"
#include "field_variable/unstructured/node_to_dof_mapping.h"
#include "basis_on_mesh/05_basis_on_mesh.h"
#include "basis_on_mesh/04_basis_on_mesh_nodes.h"
#include "mesh/unstructured_deformable.h"
#include "field_variable/field_variable_set_get.h"

namespace FieldVariable
{

/** FieldVariable class for UnstructuredDeformable mesh
 */
template<int D, typename BasisFunctionType>
class FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> :
  public FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>
{ 
public:
  //! inherited constructor 
  using FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::FieldVariableData;
  
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
  
  //! for a specific component, get all values  
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValues(std::string component, std::vector<double> &values, bool onlyNodalValues=false);
  
  //! for a specific component, get values from their global dof no.s
  template<int N>
  void getValues(std::string component, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values);
  
  //! get values from their global dof no.s for all components
  template<int N, int nComponents>
  void getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
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
  }
  
  //! for a specific component, get the values corresponding to all element-local dofs
  template<int N>
  void getElementValues(std::string component, element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
  {
    this->component_[component].getElementValues(elementNo, values);
  }
  
  //! get the values corresponding to all element-local dofs for all components
  template<int N, int nComponents>
  void getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
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
  
  //! for a specific component, get a single value from global dof no.
  double getValue(std::string component, node_no_t dofGlobalNo);

  //! get a single value from global dof no. for all components
  template<std::size_t nComponents>
  std::array<double,nComponents> getValue(node_no_t dofGlobalNo);

  //! copy the values from another field variable of the same type
  void setValues(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> &rhs);
  
  //! set values for all components for dofs, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  template<std::size_t nComponents>
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES)
  {
    std::array<int,nComponents> indices;
    
    // loop over dof numbers
    int i=0;
    for (std::vector<dof_no_t>::iterator iter = dofGlobalNos.begin(); iter != dofGlobalNos.end(); iter++, i++)
    {  
      dof_no_t dofGlobalNo = *iter;
      
      // prepare lookup indices for PETSc vector values_
      for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
      {
        indices[componentIndex] = dofGlobalNo*this->nComponents_ + componentIndex;
      }
      
      VecSetValues(this->values_, nComponents, indices.data(), values[i].data(), petscInsertMode);
    }
    
    // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called 
  }

  //! set values for dofs with a single component, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! set a single dof (all components) , after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  template<std::size_t nComponents>
  void setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value, InsertMode petscInsertMode=INSERT_VALUES)
  {
    std::array<int,nComponents> indices;
    
    // prepare lookup indices for PETSc vector values_
    for (int componentIndex = 0; componentIndex < nComponents; componentIndex++)
    {
      indices[componentIndex] = dofGlobalNo*this->nComponents_ + componentIndex;
    }
    
    VecSetValues(this->values_, nComponents, indices.data(), value.data(), petscInsertMode);
    // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called 
  }
  
  //! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
  void flushSetValues();
};

};  // namespace

#include "field_variable/unstructured/02_field_variable_set_get_unstructured_deformable.tpp"

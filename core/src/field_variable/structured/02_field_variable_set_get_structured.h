#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>

#include "field_variable/structured/01_field_variable_data_structured.h"

namespace FieldVariable
{

/** Field variable for a structured mesh, i.e. dof and node information are purely implicit.
 *  This is used for RegularFixed and StructuredDeformable meshes.
 */
template<typename BasisOnMeshType>
class FieldVariableSetGetStructured : 
  public FieldVariableDataStructured<BasisOnMeshType>
{
public:
  //! inherited constructor 
  using FieldVariableDataStructured<BasisOnMeshType>::FieldVariableDataStructured;
 
  //! for a specific component, get all values
  void getValues(std::string component, std::vector<double> &values);
  
  //! for a specific component, get values from their global dof no.s
  template<int N>
  void getValues(std::string component, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values);
  
  //! get values from their global dof no.s for all components
  template<int N, int nComponents>
  void getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
  {
    std::array<int,N*nComponents> indices;
    std::array<double,N*nComponents> result;
    
    // prepare lookup indices for PETSc vector values_
    int j=0;
    for (int i=0; i<N; i++)
    {
      for (int componentIndex = 0; componentIndex < this->nComponents_; componentIndex++, j++)
      {
        indices[j] = dofGlobalNo[i]*nComponents + componentIndex;
      }
    }
    
    VecGetValues(this->values_, N*nComponents, indices.data(), result.data());
    
    // copy result to output values
    for (int i=0; i<N; i++)
    {
      for (int componentIndex = 0; componentIndex < this->nComponents_; componentIndex++, j++)
      {
        values[i][componentIndex] = result[i*nComponents+componentIndex];
      }
    }
  }
  //! for a specific component, get the values corresponding to all element-local dofs
  template<int N>
  void getElementValues(std::string component, element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values);
  
  //! get the values corresponding to all element-local dofs for all components
  template<std::size_t nComponents>
  void getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
  {
    const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
    std::array<int,nDofsPerElement*nComponents> indices;
    std::array<double,nDofsPerElement*nComponents> result;
    
    VLOG(2) << "getElementValues element " << elementNo << ", nComponents=" << nComponents << ", " << this->nComponents_;
    
    // prepare lookup indices for PETSc vector values_
    int j=0;
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      for (int componentIndex = 0; componentIndex < this->nComponents_; componentIndex++, j++)
      {
        indices[j] = this->mesh_->getDofNo(elementNo,dofIndex)*nComponents + componentIndex;
      }
    }
    
    VecGetValues(this->values_, nDofsPerElement*nComponents, indices.data(), result.data());
    
    // copy result to output values
    for (int dofIndex=0; dofIndex<nDofsPerElement; dofIndex++)
    {
      for (int componentIndex = 0; componentIndex < this->nComponents_; componentIndex++, j++)
      {
        values[dofIndex][componentIndex] = result[dofIndex*nComponents+componentIndex];
        VLOG(2) << "getElementValues element " << elementNo << ", dofIndex " << dofIndex << " componentIndex " << componentIndex << " value: " << values[dofIndex][componentIndex];
      }
    }
  }
  //! for a specific component, get a single value from global dof no.
  double getValue(std::string component, node_no_t dofGlobalNo);

  //! get a single value from global dof no. for all components
  template<std::size_t nComponents>
  std::array<double,nComponents> getValue(node_no_t dofGlobalNo);

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

#include "field_variable/structured/02_field_variable_set_get_structured.tpp"

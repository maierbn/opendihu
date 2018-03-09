#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "data_management/finite_elements.h"
#include "output_writer/python/python_base.h"

namespace OutputWriter
{

 /*
class Python : public Generic
{
public:
 
  //! constructor
  Python(PyObject *specificSettings);
 
  //! write out solution to given filename, if timeStepNo is not -1, this value will be part of the filename
  template<typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1);
 
private:
 
  std::string filenameBase_;
};
*/


template<typename BasisOnMeshType>
class Python
{};

/* Specialization for RegularFixed. This also outputs rhs matrix and stiffness matrix for laplace problems.
 * 
 */
template<int D, typename BasisFunctionType>
class Python<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>> :
  public PythonBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> BasisOnMeshType;
 
  //! call python callback
  static PyObject *buildPyDataObject(std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables, int timeStepNo, double currentTime);  
};

// specialization for StructuredDeformable
template<int D, typename BasisFunctionType>
class Python<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> :
  public PythonBase<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
 
  //! call python callback
  static PyObject *buildPyDataObject(std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables, int timeStepNo, double currentTime);  
};

// specialization for UnstructuredDeformable
template<int D, typename BasisFunctionType>
class Python<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> :
  public PythonBase<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
 
  //! call python callback
  static PyObject *buildPyDataObject(std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables, int timeStepNo, double currentTime);  
};

};  // namespace

#include "output_writer/python/python_structured_regular_fixed.tpp"
#include "output_writer/python/python_structured_deformable.tpp"
#include "output_writer/python/python_unstructured_deformable.tpp"

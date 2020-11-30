#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "data_management/finite_element_method/finite_elements.h"
#include "output_writer/python/python_base.h"

namespace OutputWriter
{

/** Helper class that creates a python object out of a tuple of field variables.
 *  FieldVariablesForOutputWriterType is a std::tuple<std::shared_ptr<>, std::shared_ptr<>, ...> of field variables that will be output.
  */
template<typename FunctionSpaceType, typename FieldVariablesForOutputWriterType>
class Python
{};

/* Specialization for RegularFixed. This also outputs rhs matrix and stiffness matrix for laplace problems.
 *
 */
template<int D, typename BasisFunctionType, typename FieldVariablesForOutputWriterType>
class Python<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,FieldVariablesForOutputWriterType> :
  public PythonBase<FieldVariablesForOutputWriterType>
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! call python callback
  static PyObject *buildPyDataObject(FieldVariablesForOutputWriterType fieldVariables,
                                     std::string meshName, int timeStepNo, double currentTime, bool onlyNodalValues);
};

// specialization for StructuredDeformable
template<int D, typename BasisFunctionType, typename FieldVariablesForOutputWriterType>
class Python<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,FieldVariablesForOutputWriterType> :
  public PythonBase<FieldVariablesForOutputWriterType>
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! call python callback
  static PyObject *buildPyDataObject(FieldVariablesForOutputWriterType fieldVariables,
                                     std::string meshName, int timeStepNo, double currentTime, bool onlyNodalValues);
};

// specialization for Composite
template<int D, typename BasisFunctionType, typename FieldVariablesForOutputWriterType>
class Python<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,FieldVariablesForOutputWriterType> :
  public PythonBase<FieldVariablesForOutputWriterType>
{
public:
  // the function space is StructuredDeformable whereas the FunctionSpace of the class is CompositeOfDimension
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! call python callback
  static PyObject *buildPyDataObject(FieldVariablesForOutputWriterType fieldVariables,
                                     std::string meshName, int timeStepNo, double currentTime, bool onlyNodalValues);
};

// specialization for UnstructuredDeformable
template<int D, typename BasisFunctionType, typename FieldVariablesForOutputWriterType>
class Python<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FieldVariablesForOutputWriterType> :
  public PythonBase<FieldVariablesForOutputWriterType>
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  //! call python callback
  static PyObject *buildPyDataObject(FieldVariablesForOutputWriterType fieldVariables,
                                     std::string meshName, int timeStepNo, double currentTime, bool onlyNodalValues);
private:
  
  //! create a list of list where for each element the dofs are listed (if !onlyNodalValues) or the node numbers (if onlyNodalValues)
  static PyObject *buildPyElementalDofsObject(std::shared_ptr<Mesh::Mesh> meshBase, bool onlyNodalValues);
};

} // namespace

#include "output_writer/python/python_structured_regular_fixed.tpp"
#include "output_writer/python/python_structured_deformable.tpp"
#include "output_writer/python/python_composite.tpp"
#include "output_writer/python/python_unstructured_deformable.tpp"

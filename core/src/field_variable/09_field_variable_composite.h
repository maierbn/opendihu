#pragma once

#include "field_variable/08_field_variable_vector.h"

namespace FieldVariable
{

/** General field variable
 */
template<typename FunctionSpaceType,int nComponents>
class FieldVariableComposite :
  public FieldVariableVector<FunctionSpaceType,nComponents>
{
public:
  //! inherited constructor
  using FieldVariableVector<FunctionSpaceType,nComponents>::FieldVariableVector;
};

/** Partial specialization for field variable with composite mesh
 */
template<int D,typename BasisFunctionType,int nComponents>
class FieldVariableComposite<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableVector<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents>
{
public:

  typedef FunctionSpace::FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType> FunctionSpaceType;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType> SubFunctionSpaceType;

  //! inherited constructor
  using FieldVariableVector<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableVector;

  //! split this field variables into field variables for the sub function space and set the corresponding values
  void getSubFieldVariables(std::vector<std::shared_ptr<FieldVariable<SubFunctionSpaceType, nComponents>>> &subFieldVariables);

  //! get the sub field variable no i
  std::shared_ptr<FieldVariable<SubFunctionSpaceType,nComponents>> subFieldVariable(int i);

protected:

  //! initialize the subFieldVariables_ vector
  void initializeSubFieldVariables();

  std::vector<std::shared_ptr<FieldVariable<SubFunctionSpaceType,nComponents>>> subFieldVariables_;     //< field variables that have the values of this field variable on the submeshes
};



} // namespace

#include "field_variable/09_field_variable_composite.tpp"

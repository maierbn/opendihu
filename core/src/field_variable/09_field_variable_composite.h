#pragma once

#include "field_variable/08_field_variable_gradient.h"

namespace FieldVariable
{

/** General field variable
 */
template<typename FunctionSpaceType,int nComponents>
class FieldVariableComposite :
  public FieldVariableGradient<FunctionSpaceType,nComponents>
{
public:
  //! inherited constructor
  using FieldVariableGradient<FunctionSpaceType,nComponents>::FieldVariableGradient;
};

/** Partial specialization for field variable with composite mesh
 */
template<int D,typename BasisFunctionType,int nComponents>
class FieldVariableComposite<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableGradient<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents>
{
public:

  typedef FunctionSpace::FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType> FunctionSpaceType;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType> SubFunctionSpaceType;

  //! inherited constructor
  using FieldVariableGradient<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableGradient;

  //! split this field variables into field variables for the sub function space and set the corresponding values
  void getSubFieldVariables(std::vector<std::shared_ptr<FieldVariable<SubFunctionSpaceType, nComponents>>> &subFieldVariables);

  //! get the sub field variable no i with the recent values, this needs to be called collectively by all ranks
  std::shared_ptr<FieldVariable<SubFunctionSpaceType,nComponents>> subFieldVariable(int i);

  //! get the sub field variable no i with the last values, this does not need to be called collectively by all ranks
  std::shared_ptr<FieldVariable<SubFunctionSpaceType,nComponents>> subFieldVariableWithoutUpdate(int i);

  //! initialize the subFieldVariables_ vector or update the value if it was already initialized, such that subFieldVariable() gets the most recent values
  void updateSubFieldVariables();

  //! fill the gradient field with the gradient values in world coordinates of this field variable. This is only possible for scalar fields.
  void computeGradientField(std::shared_ptr<FieldVariable<FunctionSpaceType, FunctionSpaceType::dim()>> gradientField,
                            std::shared_ptr<FieldVariable<FunctionSpaceType,1>> jacobianConditionNumber = nullptr);

protected:

  std::vector<std::shared_ptr<FieldVariable<SubFunctionSpaceType,nComponents>>> subFieldVariables_;     //< field variables that have the values of this field variable on the submeshes
  std::vector<VecD<nComponents>> subFieldVariableValues_;     //< temporary buffer
};



} // namespace

#include "field_variable/09_field_variable_composite.tpp"

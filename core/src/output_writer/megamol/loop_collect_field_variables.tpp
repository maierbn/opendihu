#include "output_writer/megamol/loop_collect_field_variables.h"

#include <cstdlib>
#include "field_variable/field_variable.h"

namespace OutputWriter
{

namespace MegaMolLoopOverTuple
{
 
/** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, typename FunctionSpaceType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopCollectFieldVariables(const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName,
                          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                          std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  // call what to do in the loop body
  if (collectFieldVariables<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type,FunctionSpaceType>(std::get<i>(fieldVariables), meshName, geometryField, scalarFieldVariables))
    return;
  
  // advance iteration to next tuple element
  loopCollectFieldVariables<FieldVariablesForOutputWriterType, FunctionSpaceType, i+1>(fieldVariables, meshName, geometryField, scalarFieldVariables);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FunctionSpaceType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value
                        && CurrentFieldVariableType::element_type::nComponents() == 1 
                        && std::is_same<typename CurrentFieldVariableType::element_type::FunctionSpace, FunctionSpaceType>::value, bool>::type
collectFieldVariables(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                      std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  // if mesh name is not the specified meshName step over this field variable but do not exit the loop over field variables
  if (currentFieldVariable->functionSpace()->meshName() != meshName)
  {
    return false;  // do not break iteration
  }

  geometryField = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(currentFieldVariable->functionSpace()->geometryField());

  // number of components of currentFieldVariable is 1
  scalarFieldVariables.push_back(currentFieldVariable);

  return false;  // do not break iteration
}

// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FunctionSpaceType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value
                        && CurrentFieldVariableType::element_type::nComponents() != 1
                        && std::is_same<typename CurrentFieldVariableType::element_type::FunctionSpace, FunctionSpaceType>::value, bool>::type
collectFieldVariables(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                      std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  // if mesh name is not the specified meshName step over this field variable but do not exit the loop over field variables
  if (currentFieldVariable->functionSpace()->meshName() != meshName)
  {
    return false;  // do not break iteration
  }

  geometryField = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(currentFieldVariable->functionSpace()->geometryField());

  return false;  // do not break iteration
}

// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FunctionSpaceType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value
                        && !std::is_same<typename CurrentFieldVariableType::element_type::FunctionSpace, FunctionSpaceType>::value, bool>::type
collectFieldVariables(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                      std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  // FunctionSpaceType is different from the function space that is used in the currentFieldVariable, therefore, the geometry field cannot be extracted here
  return false;  // do not break iteration
}

// element i is of tuple type
template<typename TupleType, typename FunctionSpaceType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
collectFieldVariables(TupleType currentFieldVariableTuple, std::string meshName,
                      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                      std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  // call for tuple element
  loopCollectFieldVariables<TupleType,FunctionSpaceType>(currentFieldVariableTuple, meshName, geometryField, scalarFieldVariables);
  
  return false;  // do not break iteration
}

// element i is of vector type
template<typename VectorType, typename FunctionSpaceType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
collectFieldVariables(VectorType currentFieldVariableGradient, std::string meshName,
                      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                      std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  for (auto& currentFieldVariable : currentFieldVariableGradient)
  {
    // call function on all vector entries
    if (collectFieldVariables<typename VectorType::value_type>(
     currentFieldVariable, meshName, geometryField, scalarFieldVariables))
      return true; // break iteration
  }
  return false;  // do not break iteration
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template<typename CurrentFieldVariableType, typename FunctionSpaceType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
collectFieldVariables(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                      std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  const int D = CurrentFieldVariableType::element_type::FunctionSpace::dim();
  typedef typename CurrentFieldVariableType::element_type::FunctionSpace::BasisFunction BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType> SubFunctionSpaceType;
  const int nComponents = CurrentFieldVariableType::element_type::nComponents();

  typedef FieldVariable::FieldVariable<SubFunctionSpaceType, nComponents> SubFieldVariableType;

  std::vector<std::shared_ptr<SubFieldVariableType>> subFieldVariables;
  currentFieldVariable->getSubFieldVariables(subFieldVariables);

  for (auto& currentSubFieldVariable : subFieldVariables)
  {
    // call function on all vector entries
    if (collectFieldVariables<std::shared_ptr<SubFieldVariableType>,FunctionSpaceType>(currentSubFieldVariable, meshName, geometryField, scalarFieldVariables))
      return true;
  }

  return false;  // do not break iteration
}
}  // namespace MegaMolLoopOverTuple
}  // namespace OutputWriter

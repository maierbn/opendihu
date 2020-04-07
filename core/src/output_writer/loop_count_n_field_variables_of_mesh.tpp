#include "output_writer/loop_count_n_field_variables_of_mesh.h"

#include <cstdlib>
#include "field_variable/field_variable.h"

namespace OutputWriter
{

namespace LoopOverTuple
{
 
/** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopCountNFieldVariablesOfMesh(const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName,
                               int &nFieldVariablesOfMesh)
{
  // call what to do in the loop body
  if (countNFieldVariablesOfMesh<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type>(std::get<i>(fieldVariables), meshName, nFieldVariablesOfMesh))
    return;
  
  // advance iteration to next tuple element
  loopCountNFieldVariablesOfMesh<FieldVariablesForOutputWriterType, i+1>(fieldVariables, meshName, nFieldVariablesOfMesh);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
countNFieldVariablesOfMesh(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                           int &nFieldVariablesOfMesh)
{
  VLOG(2) << "count number of field variables in type " << StringUtility::demangle(typeid(CurrentFieldVariableType).name())
    << ", before: " << nFieldVariablesOfMesh;

  // if the field variable is a null pointer, return but do not break iteration
  if (!currentFieldVariable)
    return false;

  // if mesh name is the specified meshName, count this mesh
  if (currentFieldVariable->functionSpace()->meshName() == meshName)
  {
    nFieldVariablesOfMesh++;
  }
  
  return false;  // do not break iteration 
}

// element i is of vector type
template<typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
countNFieldVariablesOfMesh(VectorType currentFIeldVariableGradient, std::string meshName,
                           int &nFieldVariablesOfMesh)
{
  VLOG(2) << "count number of field variables in vector with size " << currentFIeldVariableGradient.size();
  for (auto& currentFieldVariable : currentFIeldVariableGradient)
  {
    // call function on all vector entries
    if (countNFieldVariablesOfMesh<typename VectorType::value_type>(currentFieldVariable, meshName, nFieldVariablesOfMesh))
      return true;
  }
  
  return false;  // do not break iteration
}

// element i is of tuple type
template<typename TupleType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
countNFieldVariablesOfMesh(TupleType currentFieldVariableTuple, std::string meshName,
                           int &nFieldVariablesOfMesh)
{
  // call for tuple element
  loopCountNFieldVariablesOfMesh<TupleType>(currentFieldVariableTuple, meshName, nFieldVariablesOfMesh);
 
  return false;  // do not break iteration
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template<typename CurrentFieldVariableType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
countNFieldVariablesOfMesh(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                           int &nFieldVariablesOfMesh)
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
    if (countNFieldVariablesOfMesh<std::shared_ptr<SubFieldVariableType>>(currentSubFieldVariable, meshName, nFieldVariablesOfMesh))
      return true;
  }

  return false;  // do not break iteration
}
}  // namespace LoopOverTuple
}  // namespace OutputWriter

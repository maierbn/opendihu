#include "output_writer/exfile/loop_output_exnode.h"

#include "output_writer/loop_count_n_field_variables_of_mesh.h"

#include <cstdlib>
#include "field_variable/field_variable.h"

namespace OutputWriter
{

namespace ExfileLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, typename AllFieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopOutputExnode(const FieldVariablesForOutputWriterType &fieldVariables, const AllFieldVariablesForOutputWriterType &allFieldVariables, std::string meshName, 
                 std::ofstream &file
)
{
  // call what to do in the loop body
  if (outputExnode<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type, AllFieldVariablesForOutputWriterType>(
        std::get<i>(fieldVariables), allFieldVariables, meshName, file))
    return;
  
  // advance iteration to next tuple element
  loopOutputExnode<FieldVariablesForOutputWriterType, AllFieldVariablesForOutputWriterType, i+1>(fieldVariables, allFieldVariables, meshName, file);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
outputExnode(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
             std::ofstream &file)
{
  // if mesh name is the specified meshName
  if (currentFieldVariable->functionSpace()->meshName() == meshName)
  {
    // here we have the type of the mesh with meshName (which is typedef to FunctionSpace)
    typedef typename CurrentFieldVariableType::element_type::FunctionSpace FunctionSpace;
   
    // count number of field variables for the particular mesh
    int nFieldVariablesInMesh = 0;
    LoopOverTuple::loopCountNFieldVariablesOfMesh(fieldVariables, meshName, nFieldVariablesInMesh);
    
    // call exfile writer to output all field variables with the meshName
    ExfileWriter<FunctionSpace, FieldVariablesForOutputWriterType>::outputExnode(file, fieldVariables, meshName, currentFieldVariable->functionSpace(), nFieldVariablesInMesh);
   
    return true;  // break iteration
  }
  
  return false;  // do not break iteration 
}

// Elementent i is of vector type
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
outputExnode(VectorType currentFieldVariableVector, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
             std::ofstream &file)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (outputExnode<typename VectorType::value_type,FieldVariablesForOutputWriterType>(currentFieldVariable, fieldVariables, meshName, file))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// Elementent i is of tuple type
template<typename TupleType, typename AllFieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
outputExnode(TupleType currentFieldVariableTuple, const AllFieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
             std::ofstream &file)
{
  // call for tuple element
  loopOutputExnode<TupleType, AllFieldVariablesForOutputWriterType>(currentFieldVariableTuple, fieldVariables, meshName, file);
  
  return false;  // do not break iteration 
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template<typename CurrentFieldVariableType, typename AllFieldVariablesForOutputWriterType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
outputExnode(CurrentFieldVariableType currentFieldVariable, const AllFieldVariablesForOutputWriterType &fieldVariables, std::string meshName,
             std::ofstream &file)
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
    if (outputExnode<std::shared_ptr<SubFieldVariableType>,AllFieldVariablesForOutputWriterType>(currentSubFieldVariable, fieldVariables, meshName, file))
      return true;
  }

  return false;  // do not break iteration
}
}  // namespace ExfileLoopOverTuple
}  // namespace OutputWriter


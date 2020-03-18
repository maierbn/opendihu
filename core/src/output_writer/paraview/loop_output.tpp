#include "output_writer/paraview/loop_output.h"

#include "output_writer/paraview/paraview_writer.h"
#include "output_writer/loop_count_n_field_variables_of_mesh.h"

#include <cstdlib>
#include "field_variable/field_variable.h"

namespace OutputWriter
{

namespace ParaviewLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, typename AllFieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopOutput(const FieldVariablesForOutputWriterType &fieldVariables, const AllFieldVariablesForOutputWriterType &allFieldVariables,
           std::string meshName, std::string filename, PythonConfig specificSettings
)
{
  // call what to do in the loop body
  if (output<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type, AllFieldVariablesForOutputWriterType>(
        std::get<i>(fieldVariables), allFieldVariables, meshName, filename, specificSettings))
    return;
  
  // advance iteration to next tuple element
  loopOutput<FieldVariablesForOutputWriterType, AllFieldVariablesForOutputWriterType, i+1>(fieldVariables, allFieldVariables, meshName, filename, specificSettings);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
output(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
       std::string filename, PythonConfig specificSettings)
{
  // if mesh name is the specified meshName
  if (currentFieldVariable->functionSpace()->meshName() == meshName)
  {
    // here we have the type of the mesh with meshName (which is typedef to FunctionSpace)
    typedef typename CurrentFieldVariableType::element_type::FunctionSpace FunctionSpace;
   
    int nFieldVariablesInMesh = 0;
    LoopOverTuple::loopCountNFieldVariablesOfMesh(fieldVariables, meshName, nFieldVariablesInMesh);
    
    // call exfile writer to output all field variables with the meshName
    ParaviewWriter<FunctionSpace, FieldVariablesForOutputWriterType>::outputFile(filename, fieldVariables, meshName, currentFieldVariable->functionSpace(), nFieldVariablesInMesh, specificSettings);
   
    return true;  // break iteration
  }
  
  return false;  // do not break iteration 
}

// element i is of vector type
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
output(VectorType currentFieldVariableVector, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
       std::string filename, PythonConfig specificSettings)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (output<typename VectorType::value_type,FieldVariablesForOutputWriterType>(currentFieldVariable, fieldVariables, meshName, filename, specificSettings))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename AllFieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
output(TupleType currentFieldVariableTuple, const AllFieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
       std::string filename, PythonConfig specificSettings)
{
  // call for tuple element
  loopOutput<TupleType, AllFieldVariablesForOutputWriterType>(currentFieldVariableTuple, fieldVariables, meshName, filename, specificSettings);
  
  return false;  // do not break iteration 
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template<typename CurrentFieldVariableType, typename AllFieldVariablesForOutputWriterType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
output(CurrentFieldVariableType currentFieldVariable, const AllFieldVariablesForOutputWriterType &fieldVariables, std::string meshName,
       std::string filename, PythonConfig specificSettings)
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
    if (output<std::shared_ptr<SubFieldVariableType>,AllFieldVariablesForOutputWriterType>(currentSubFieldVariable, fieldVariables, meshName, filename, specificSettings))
      return true;
  }

  return false;  // do not break iteration
}
}  // namespace ParaviewLoopOverTuple
}  // namespace OutputWriter

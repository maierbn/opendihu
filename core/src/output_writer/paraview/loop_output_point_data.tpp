#include "output_writer/paraview/loop_output_point_data.h"

#include <cstdlib>

#include "field_variable/field_variable.h"
#include "output_writer/paraview/paraview.h"

namespace OutputWriter
{

// forward declaration
struct Paraview;
 
namespace ParaviewLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopOutputPointData(const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
                    std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement
)
{
  // call what to do in the loop body
  if (outputPointData<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type, FieldVariablesForOutputWriterType>(
        std::get<i>(fieldVariables), fieldVariables, meshName, file, binaryOutput, fixedFormat, onlyParallelDatasetElement))
    return;
  
  // advance iteration to next tuple element
  loopOutputPointData<FieldVariablesForOutputWriterType, i+1>(fieldVariables, meshName, file, binaryOutput, fixedFormat, onlyParallelDatasetElement);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
outputPointData(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
                std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement)
{
  // if mesh name is the specified meshName
  if (currentFieldVariable->functionSpace()->meshName() == meshName && !currentFieldVariable->isGeometryField())
  {
    Paraview::writeParaviewFieldVariable<typename CurrentFieldVariableType::element_type>(*currentFieldVariable, file, binaryOutput, fixedFormat, onlyParallelDatasetElement);
  }
  
  return false;  // do not break iteration 
}

// Elementent i is of vector type
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
outputPointData(VectorType currentFIeldVariableGradient, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
                std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement)
{
  for (auto& currentFieldVariable : currentFIeldVariableGradient)
  {
    // call function on all vector entries
    if (outputPointData<typename VectorType::value_type,FieldVariablesForOutputWriterType>(currentFieldVariable, fieldVariables, meshName, file, binaryOutput, fixedFormat, onlyParallelDatasetElement))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
outputPointData(TupleType currentFieldVariableTuple, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
                std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement)
{
  // call for tuple element
  loopOutputPointData<TupleType>(currentFieldVariableTuple, meshName, file, binaryOutput, fixedFormat, onlyParallelDatasetElement);
  
  return false;  // do not break iteration 
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
outputPointData(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName,
                std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement)
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
    if (outputPointData<std::shared_ptr<SubFieldVariableType>,FieldVariablesForOutputWriterType>(currentSubFieldVariable, fieldVariables, meshName, file, binaryOutput, fixedFormat, onlyParallelDatasetElement))
      return true;
  }

  return false;  // do not break iteration
}
}  // namespace ParaviewLoopOverTuple
}  // namespace OutputWriter

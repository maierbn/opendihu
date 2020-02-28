#include "output_writer/paraview/loop_collect_field_variables_names.h"

#include "output_writer/paraview/paraview_writer.h"
#include "easylogging++.h"
#include "utility/vector_operators.h"
#include <cstdlib>
#include "field_variable/field_variable.h"

namespace OutputWriter
{

namespace ParaviewLoopOverTuple
{

 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopGetGeometryFieldNodalValues(const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> meshNames,
                   std::vector<double> &values
)
{
  // call what to do in the loop body
  if (getGeometryFieldNodalValues<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type, FieldVariablesForOutputWriterType>(
        std::get<i>(fieldVariables), fieldVariables, meshNames, values))
    return;

  // advance iteration to next tuple element
  loopGetGeometryFieldNodalValues<FieldVariablesForOutputWriterType, i+1>(fieldVariables, meshNames, values);
}

// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
getGeometryFieldNodalValues(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> meshNames,
               std::vector<double> &values)
{
  VLOG(1) << "getGeometryFieldNodalValues meshNames: " << meshNames
    << ", own: " << currentFieldVariable->functionSpace()->meshName()
    << ", fieldVariable name \"" << currentFieldVariable->name() << "\""
    << ", isGeometryField: " << currentFieldVariable->isGeometryField();

  // if mesh name is one of the specified meshNames and it is the geometry field
  if (meshNames.find(currentFieldVariable->functionSpace()->meshName()) != meshNames.end()
    && currentFieldVariable->isGeometryField())
  {
    const int nComponents = CurrentFieldVariableType::element_type::nComponents();
    std::array<std::vector<double>, nComponents> componentValues;

    // initialize the dofNosLocalNaturalOrdering vector of the meshPartition to be able to get the values in the natural ordering
    currentFieldVariable->functionSpace()->meshPartition()->initializeDofNosLocalNaturalOrdering(currentFieldVariable->functionSpace());

    VLOG(1) << "nComponents: " << nComponents;

    // get all local values with ghosts for the components
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      std::vector<double> retrievedLocalValues;
      currentFieldVariable->getValues(componentNo, currentFieldVariable->functionSpace()->meshPartition()->dofNosLocalNaturalOrdering(),
                              retrievedLocalValues);

      const int nDofsPerNode = CurrentFieldVariableType::element_type::FunctionSpace::nDofsPerNode();
      const node_no_t nNodesLocal = currentFieldVariable->functionSpace()->meshPartition()->nNodesLocalWithGhosts();

      VLOG(1) << "nNodesLocal: " << nNodesLocal;

      // for Hermite only extract the non-derivative values
      componentValues[componentNo].resize(nNodesLocal);

      int index = 0;
      for (int i = 0; i < nNodesLocal; i++)
      {
        componentValues[componentNo][i] = retrievedLocalValues[index];
        index += nDofsPerNode;
      }
    }
    values.reserve(values.size()+componentValues[0].size()*nComponents);

    // copy values in consecutive order (x y z x y z) to output
    for (int i = 0; i < componentValues[0].size(); i++)
    {
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        values.push_back(componentValues[componentNo][i]);
      }
    }
  }

  return false;  // do not break iteration
}

// element i is of vector type
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
getGeometryFieldNodalValues(VectorType currentFieldVariableVector, const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> meshNames,
               std::vector<double> &values)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (getGeometryFieldNodalValues<typename VectorType::value_type,FieldVariablesForOutputWriterType>(currentFieldVariable, fieldVariables, meshNames, values))
      return true; // break iteration
  }
  return false;  // do not break iteration
}

// element i is of tuple type
template<typename TupleType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
getGeometryFieldNodalValues(TupleType currentFieldVariableTuple, const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> meshNames,
                            std::vector<double> &values)
{
  // call for tuple element
  loopGetGeometryFieldNodalValues<TupleType>(currentFieldVariableTuple, meshNames, values);

  return false;  // do not break iteration
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
getGeometryFieldNodalValues(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> meshNames,
                            std::vector<double> &values)
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
    if (getGeometryFieldNodalValues<std::shared_ptr<SubFieldVariableType>,FieldVariablesForOutputWriterType>(currentSubFieldVariable, fieldVariables, meshNames, values))
      return true;
  }

  return false;  // do not break iteration
}
}  // namespace ParaviewLoopOverTuple
}  // namespace OutputWriter

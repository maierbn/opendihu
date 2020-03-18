#include "output_writer/paraview/loop_collect_field_variables_names.h"

#include "output_writer/paraview/paraview_writer.h"

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
loopGetNodalValues(const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> meshNames,
                   std::map<std::string,std::vector<double>> &values
)
{
  // call what to do in the loop body
  if (getNodalValues<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type, FieldVariablesForOutputWriterType>(
        std::get<i>(fieldVariables), fieldVariables, meshNames, values))
    return;
  
  // advance iteration to next tuple element
  loopGetNodalValues<FieldVariablesForOutputWriterType, i+1>(fieldVariables, meshNames, values);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
getNodalValues(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> meshNames,
               std::map<std::string,std::vector<double>> &values)
{
  VLOG(1) << "field variable " << StringUtility::demangle(typeid(currentFieldVariable).name()) << " name \"" << currentFieldVariable->name()
    << "\", is geometry: " << currentFieldVariable->isGeometryField() << ", values size: " << values.size();

  if (!currentFieldVariable->functionSpace())
    return false;

  // if mesh name is one of the specified meshNames (and it is not a geometry field)
  if (meshNames.find(currentFieldVariable->functionSpace()->meshName()) != meshNames.end()
    && !currentFieldVariable->isGeometryField())
  {
    const int nComponents = CurrentFieldVariableType::element_type::nComponents();
    std::array<std::vector<double>, nComponents> componentValues;

    // initialize the dofNosLocalNaturalOrdering vector of the meshPartition to be able to get the values in the natural ordering
    currentFieldVariable->functionSpace()->meshPartition()->initializeDofNosLocalNaturalOrdering(currentFieldVariable->functionSpace());

    // ensure that ghost values are in place
    currentFieldVariable->zeroGhostBuffer();
    currentFieldVariable->setRepresentationGlobal();
    currentFieldVariable->startGhostManipulation();

    // get all local values including ghosts for the components
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      std::vector<double> retrievedLocalValues;
      const std::vector<dof_no_t> &dofNosLocalNaturalOrdering = currentFieldVariable->functionSpace()->meshPartition()->dofNosLocalNaturalOrdering();
      currentFieldVariable->getValues(componentNo, dofNosLocalNaturalOrdering, retrievedLocalValues);

      const int nDofsPerNode = CurrentFieldVariableType::element_type::FunctionSpace::nDofsPerNode();
      const node_no_t nNodesLocal = currentFieldVariable->functionSpace()->meshPartition()->nNodesLocalWithGhosts();

      // for Hermite only extract the non-derivative values
      componentValues[componentNo].resize(nNodesLocal);

      int index = 0;
      for (int i = 0; i < nNodesLocal; i++)
      {
        componentValues[componentNo][i] = retrievedLocalValues[index];
        index += nDofsPerNode;
      }
    }

    std::string fieldVariableName = currentFieldVariable->name();

    // create entry for field variable name if it does not exist and reserve enough space for all values
    values[fieldVariableName].reserve(values[fieldVariableName].size() + componentValues[0].size()*nComponents);
    LOG(DEBUG) << "get nodal values of \"" << fieldVariableName << "\" (" << currentFieldVariable << ").";

    // copy values in consecutive order (x y z x y z) to output
    for (int i = 0; i < componentValues[0].size(); i++)
    {
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        values[fieldVariableName].push_back(componentValues[componentNo][i]);
      }
    }

  }
  
  return false;  // do not break iteration 
}

// element i is of vector type
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
getNodalValues(VectorType currentFieldVariableVector, const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> meshNames,
               std::map<std::string,std::vector<double>> &values)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (getNodalValues<typename VectorType::value_type,FieldVariablesForOutputWriterType>(currentFieldVariable, fieldVariables, meshNames, values))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
getNodalValues(TupleType currentFieldVariableTuple, const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> meshNames,
               std::map<std::string,std::vector<double>> &values)
{
  // call for tuple element
  loopGetNodalValues<TupleType>(currentFieldVariableTuple, meshNames, values);
 
  return false;  // do not break iteration 
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
getNodalValues(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> meshNames,
               std::map<std::string,std::vector<double>> &values)
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
    if (getNodalValues<std::shared_ptr<SubFieldVariableType>,FieldVariablesForOutputWriterType>(currentSubFieldVariable, fieldVariables, meshNames, values))
      return true;
  }

  return false;  // do not break iteration
}
}  // namespace ParaviewLoopOverTuple
}  // namespace OutputWriter

#include "output_writer/exfile/loop_output_node_values.h"

#include <cstdlib>

namespace OutputWriter
{

namespace ExfileLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopOutputNodeValues(const OutputFieldVariablesType &fieldVariables, std::string meshName,
                     std::ostream &stream, node_no_t nodeGlobalNo)
{
  // call what to do in the loop body
  if (outputNodeValues<typename std::tuple_element<i,OutputFieldVariablesType>::type>(
    std::get<i>(fieldVariables), meshName, stream, nodeGlobalNo))
    return;
  
  // advance iteration to next tuple element
  loopOutputNodeValues<OutputFieldVariablesType, i+1>(fieldVariables, meshName, stream, nodeGlobalNo);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType>
typename std::enable_if<!TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
outputNodeValues(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                      std::ostream &stream, node_no_t nodeGlobalNo)
{
  VLOG(2) << "loop_output_node_values.tpp:34, outputNodeValues, field variable " << currentFieldVariable->name()<<" at " << currentFieldVariable->values();
 
  // if mesh name is not the specified meshName step over this field variable but do not exit the loop over field variables
  if (currentFieldVariable->mesh()->meshName() != meshName)
  {
    return false;  // do not break iteration
  }
  
  // output the values of a node for the current field variable

  // get the class of the current field variable
  typedef CurrentFieldVariableType FieldVariableType;

  // get the current field variable
  auto fieldVariable = currentFieldVariable;  // this type is std::shared_ptr<FieldVariable<..>>

  const int nDofsPerNode = FieldVariableType::element_type::BasisOnMesh::nDofsPerNode();

  // get all dofs that are associated with the current node
  //std::array<dof_no_t,nDofsPerNode> dofGlobalNos;
  std::vector<dof_no_t> dofGlobalNos;
  dofGlobalNos.reserve(nDofsPerNode);
  fieldVariable->mesh()->getNodeDofs(nodeGlobalNo, dofGlobalNos);

  // loop over components of the field variable
  for (int componentNo = 0; componentNo < fieldVariable->nComponents(); componentNo++)
  {
    // get the values of the dofs for the current component
    //std::array<double,nDofsPerNode> values;
    //fieldVariable->template getValues<nDofsPerNode>(componentNo, dofGlobalNos, values);

    std::vector<double> values;
    values.reserve(nDofsPerNode);
    
    VLOG(2) << "get dofGlobalNos: " << dofGlobalNos;
    VLOG(2) << "field variable is geometry field: " << fieldVariable->isGeometryField();
    
    fieldVariable->getValues(componentNo, dofGlobalNos, values);
    
    // output values
    for (double value : values)
    {
      stream << "  " << std::scientific << std::setprecision(17) << value << std::endl;
    }
  }
  
  return false;  // do not break iteration 
}

// element i is of vector type
template<typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
outputNodeValues(VectorType currentFieldVariableVector, std::string meshName,
                 std::ostream &stream, node_no_t nodeGlobalNo)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (outputNodeValues<typename VectorType::value_type>(currentFieldVariable, meshName, stream, nodeGlobalNo))
      return true;
  }
  
  return false;  // do not break iteration
}

};  //namespace ExfileLoopOverTuple
};  //namespace OutputWriter

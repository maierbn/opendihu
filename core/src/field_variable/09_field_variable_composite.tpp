#include "field_variable/09_field_variable_composite.h"

namespace FieldVariable
{

//! split this field variables into field variables for the sub function space and set the corresponding values
template<int D,typename BasisFunctionType,int nComponents>
void FieldVariableComposite<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,nComponents>::
getSubFieldVariables(std::vector<std::shared_ptr<FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>, nComponents>>> &subFieldVariables)
{
  subFieldVariables.clear();

  // get subFunctionSpaces
  std::vector<std::shared_ptr<SubFunctionSpaceType>> subFunctionSpaces = this->functionSpace_->subFunctionSpaces();

  // get own values
  std::vector<VecD<nComponents>> ownValues;
  this->getValuesWithoutGhosts(ownValues);

  // loop over sub function spaces
  int subMeshNo = 0;
  for (std::shared_ptr<SubFunctionSpaceType> &subFunctionSpace : subFunctionSpaces)
  {
    // determine name of sub field variable
    std::stringstream fieldVariableName;
    fieldVariableName << this->name_ << "_" << subFunctionSpace->meshName();

    // create sub field variable
    std::vector<std::string> componentNames(this->componentNames_.begin(), this->componentNames_.end());
    std::shared_ptr<FieldVariable<SubFunctionSpaceType,nComponents>> newFieldVariable = subFunctionSpace->template createFieldVariable<nComponents>(fieldVariableName.str(), componentNames);

    // determine values for the sub field variable
    std::vector<VecD<nComponents>> subFieldVariableValues(subFunctionSpace->nDofsLocalWithoutGhosts());

    // loop over dofs
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < subFunctionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t compositeNodeNo = this->functionSpace_->meshPartition()->getNodeNoLocalFromSubmesh(subMeshNo, nodeNoLocal);

      for (int nodalDofNo = 0; nodalDofNo < subFunctionSpace->nDofsPerNode(); nodalDofNo++)
      {
        dof_no_t compositeDofNo = compositeNodeNo*subFunctionSpace->nDofsPerNode() + nodalDofNo;
        dof_no_t subFunctionSpaceDofNo = nodeNoLocal*subFunctionSpace->nDofsPerNode() + nodalDofNo;

        subFieldVariableValues[subFunctionSpaceDofNo] = ownValues[compositeDofNo];
      }
    }

    // set values in the sub field variable
    newFieldVariable->startGhostManipulation();
    newFieldVariable->setValuesWithoutGhosts(subFieldVariableValues);
    newFieldVariable->finishGhostManipulation();

    // add field variable to subFieldVariables
    subFieldVariables.push_back(newFieldVariable);
    subMeshNo++;
  }
}

} // namespace

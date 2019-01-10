#include "output_writer/megamol/megamol_writer.h"

#include "output_writer/megamol/loop_collect_field_variables.h"

namespace OutputWriter
{

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
void MegaMOLWriter<FunctionSpaceType, OutputFieldVariablesType>::
outputData(OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<FunctionSpaceType> functionSpace,
           PythonConfig specificSettings)
{
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField;
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> scalarFieldVariables;

  // collect the geometryField and all scalar field variables for the current mesh
  MegaMOLLoopOverTuple::loopCollectFieldVariables(fieldVariables, meshName, geometryField, scalarFieldVariables);

  LOG(DEBUG) << "mesh \"" << meshName << "\", has " << scalarFieldVariables.size() << " scalar field variables";

  // retrieve the geometry field values
  std::vector<Vec3> geometryFieldValues;
  geometryField->getValuesWithoutGhosts(geometryFieldValues);

  LOG(DEBUG) << geometryFieldValues.size() << " geometryField values: " << geometryFieldValues;

  // get all other scalar field variables
  std::vector<std::vector<double>> values(scalarFieldVariables.size());
  for (int i = 0; i < scalarFieldVariables.size(); i++)
  {
    scalarFieldVariables[i]->getValuesWithoutGhosts(values[i]);
    LOG(DEBUG) << "field variable \"" << scalarFieldVariables[i]->name() << "\" has " << values[i].size() << " values: " << values[i];
  }
}
  

};

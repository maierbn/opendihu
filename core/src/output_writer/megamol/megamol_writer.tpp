#include "output_writer/megamol/megamol_writer.h"

#include "output_writer/megamol/loop_collect_field_variables.h"
#include <adios2.h>

namespace OutputWriter
{

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
std::map<std::string, adios2::Variable<double>> MegaMolWriter<FunctionSpaceType,OutputFieldVariablesType>::geometryTable_;

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
std::map<std::string, adios2::Variable<double>> MegaMolWriter<FunctionSpaceType,OutputFieldVariablesType>::vmTable_;

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
void MegaMolWriter<FunctionSpaceType, OutputFieldVariablesType>::
outputData(OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<FunctionSpaceType> functionSpace,
           PythonConfig specificSettings, std::shared_ptr<adios2::Engine> adiosWriter, std::shared_ptr<adios2::IO> adiosIo)
{
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField;
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> scalarFieldVariables;

  // collect the geometryField and all scalar field variables for the current mesh
  MegaMolLoopOverTuple::loopCollectFieldVariables(fieldVariables, meshName, geometryField, scalarFieldVariables);

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

  // determine the number of values on the local rank
  const long unsigned int nValuesGeometryTable = (long unsigned int) geometryFieldValues.size()*3;

  // create ADIOS field variable
  if (geometryTable_.find(meshName) == geometryTable_.end())
  {
    geometryTable_[meshName] = adiosIo->DefineVariable<double>(meshName, {adios2::JoinedDim, nValuesGeometryTable}, {}, {1, nValuesGeometryTable});

    if (!values.empty())
    {
      vmTable_[meshName] = adiosIo->DefineVariable<double>(meshName, {adios2::JoinedDim, nValuesGeometryTable}, {}, {1, nValuesGeometryTable});
    }
  }

  // convert data to be send to ADIOS
  // geometryTable
  std::vector<double> geometryTableData(nValuesGeometryTable);
  for (int i = 0; i < geometryFieldValues.size(); i++)
  {
    for (int j = 0; j != 3; j++)
    {
      geometryTableData[3*i + j] = geometryFieldValues[i][j];
    }
  }

  // vmTable
  if (!values.empty())
  {
    adiosWriter->Put<double>(vmTable_[meshName], values[0].data());
  }

  adiosWriter->Put<double>(geometryTable_[meshName], geometryTableData.data());

}
  
}  // namespace

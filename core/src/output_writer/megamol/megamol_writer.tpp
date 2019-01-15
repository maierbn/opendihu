#include "output_writer/megamol/megamol_writer.h"

#include "output_writer/megamol/loop_collect_field_variables.h"

#ifdef HAVE_ADIOS
#include <adios2.h>
#endif

namespace OutputWriter
{

#ifdef HAVE_ADIOS
template<typename FunctionSpaceType, typename OutputFieldVariablesType>
std::map<std::string, adios2::Variable<double>> MegaMolWriter<FunctionSpaceType,OutputFieldVariablesType>::geometryTable_;

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
std::map<std::string, adios2::Variable<double>> MegaMolWriter<FunctionSpaceType,OutputFieldVariablesType>::vmTable_;

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
void MegaMolWriter<FunctionSpaceType, OutputFieldVariablesType>::
outputData(OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<FunctionSpaceType> functionSpace,
           PythonConfig specificSettings, std::shared_ptr<adios2::Engine> adiosWriter, std::shared_ptr<adios2::IO> adiosIo, BoundingBox &boundingBox)
{
  LOG(DEBUG) << "MegaMolWriter::outputData";

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

  // create ADIOS fielddefine variable variable
  if (geometryTable_.find(meshName) == geometryTable_.end())
  {
    std::stringstream variableName;
    //variableName << meshName << "_geometry";
    variableName << "xyz";


    int localSize = nValuesGeometryTable;
    int offset = 0;
    MPI_Exscan(&localSize, &offset, 1, MPI_INT, MPI_SUM, functionSpace->meshPartition()->mpiCommunicator());
    LOG(DEBUG) << "localSize: " << localSize << ", offset: " << offset << ", globalSize: " << geometryField->nDofsGlobal()*3;

    LOG(DEBUG) << "define variable \"" << variableName.str() << "\". nValuesGeometryTable: " << nValuesGeometryTable;
    geometryTable_[meshName] = adiosIo->DefineVariable<double>(variableName.str(), {(long unsigned int)geometryField->nDofsGlobal()*3}, {(long unsigned int)offset}, {(long unsigned int)nValuesGeometryTable});

    if (!values.empty())
    {
      std::stringstream variableName;
      //variableName << meshName << "_vm";
      variableName << "i";
      LOG(DEBUG) << "define variable \"" << variableName.str() << "\". size: " << values[0].size();
      vmTable_[meshName] = adiosIo->DefineVariable<double>(variableName.str(), {adios2::JoinedDim}, {}, {values[0].size()});
    }
  }
  else
  {
    LOG(DEBUG) << "geometryTable[" << meshName << "] is already set:";
    for(std::map<std::string, adios2::Variable<double>>::iterator iter = geometryTable_.begin(); iter != geometryTable_.end(); iter++)
    {
      LOG(DEBUG) << "key " << iter->first;
    }
  }

  // compute local bounding box
  if (!boundingBox.initialized)
  {
    boundingBox.min = geometryFieldValues[0];
    boundingBox.max = geometryFieldValues[0];
    boundingBox.initialized = true;
  }
  for(std::vector<Vec3>::iterator iter = geometryFieldValues.begin(); iter != geometryFieldValues.end(); iter++)
  {
    for (int i = 0; i < 3; i++)
    {
      if ((*iter)[i] < boundingBox.min[i])
      {
        boundingBox.min[i] = (*iter)[i];
      }
      if ((*iter)[i] > boundingBox.max[i])
      {
        boundingBox.max[i] = (*iter)[i];
      }
    }
  }

  LOG(DEBUG) << "local min: " << boundingBox.min << ", local max: " << boundingBox.max;

  // convert data to be send to ADIOS
  // geometryTable
  std::vector<double> geometryTableData(nValuesGeometryTable);
  double maximum = 0;
  for (int i = 0; i < geometryFieldValues.size(); i++)
  {
    for (int j = 0; j != 3; j++)
    {
      geometryTableData[3*i + j] = geometryFieldValues[i][j];
      maximum = std::max(maximum, geometryFieldValues[i][j]);
    }
  }
  LOG(DEBUG) << "maximum of geometryFieldValues: " << maximum;

  LOG(DEBUG) << "nValuesGeometryTable: " << nValuesGeometryTable << ", geometryTableData: " << geometryTableData;
  LOG(DEBUG) << "values[0].size(): " << values[0].size() << ", values[0]: " << values[0];

  // vmTable
  if (!values.empty())
  {
    adiosWriter->Put<double>(vmTable_[meshName], values[0].data());
  }

  adiosWriter->Put<double>(geometryTable_[meshName], geometryTableData.data());

}
#endif

}  // namespace

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
std::map<std::string, std::shared_ptr<std::vector<double>>> MegaMolWriter<FunctionSpaceType,OutputFieldVariablesType>::geometryTableData_;

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
std::map<std::string, std::shared_ptr<std::vector<double>>> MegaMolWriter<FunctionSpaceType,OutputFieldVariablesType>::vmTableData_;

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
std::shared_ptr<std::vector<Vec3>> MegaMolWriter<FunctionSpaceType,OutputFieldVariablesType>::geometryFieldValues_;

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
void MegaMolWriter<FunctionSpaceType, OutputFieldVariablesType>::
outputData(OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<FunctionSpaceType> functionSpace,
           PythonConfig specificSettings, std::shared_ptr<adios2::Engine> adiosWriter, std::shared_ptr<adios2::IO> adiosIo,
           BoundingBox &boundingBox, double &minimalDistanceBetweenFibers)
{
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField;
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> scalarFieldVariables;

  // collect the geometryField and all scalar field variables for the current mesh
  MegaMolLoopOverTuple::loopCollectFieldVariables(fieldVariables, meshName, geometryField, scalarFieldVariables);

  // retrieve the new geometry field values
  std::shared_ptr<std::vector<Vec3>> newGeometryFieldValues = std::make_shared<std::vector<Vec3>>();
  geometryField->getValuesWithoutGhosts(*newGeometryFieldValues);

  // compute the minimal distance between two neighbouring fibers, by comparing all points of the current fiber with the points of the previous fiber
  // the previous fiber data is still in geometryFieldValues_
  if (geometryFieldValues_ && geometryFieldValues_ != newGeometryFieldValues)
  {
    if (!geometryFieldValues_->empty() && !newGeometryFieldValues->empty())
    {
      // initialize minimalDistanceBetweenFibers the first time
      double distanceBetweenPoints = MathUtility::distance<3>((*geometryFieldValues_)[0], (*newGeometryFieldValues)[0]);
      if (minimalDistanceBetweenFibers == -1.0)
      {
        LOG(DEBUG) << "initialize minimalDistanceBetweenFibers to " << distanceBetweenPoints;
        minimalDistanceBetweenFibers = distanceBetweenPoints;
      }
    }
    for (int i = 1; i < std::min(geometryFieldValues_->size(), newGeometryFieldValues->size()); i++)
    {
      double distanceBetweenPoints = MathUtility::distance<3>((*geometryFieldValues_)[i], (*newGeometryFieldValues)[i]);
      minimalDistanceBetweenFibers = std::min(minimalDistanceBetweenFibers, distanceBetweenPoints);
      LOG(DEBUG) << "  i: " << i << ", minimalDistanceBetweenFibers: " << minimalDistanceBetweenFibers << ", distanceBetweenPoints: " << distanceBetweenPoints;
    }
  }
  else
  {
    LOG(DEBUG) << "mesh: " << meshName << ", geometryFieldValues_ is not set or there is only one fiber per rank";
  }

  // update geometryFieldValues_ with the new values
  geometryFieldValues_ = newGeometryFieldValues;

  std::shared_ptr<Partition::RankSubset> rankSubset = DihuContext::partitionManager()->rankSubsetForCollectiveOperations();

  LOG(DEBUG) << "MegaMolWriter::outputDatam, mesh \"" << meshName << "\", has " << scalarFieldVariables.size() << " scalar field variables";
  LOG(DEBUG) << geometryFieldValues_->size() << " geometryField values: " << geometryFieldValues_;

  // get all other scalar field variables
  std::vector<std::vector<double>> values(scalarFieldVariables.size());
  for (int i = 0; i < scalarFieldVariables.size(); i++)
  {
    scalarFieldVariables[i]->getValuesWithoutGhosts(values[i]);
    LOG(DEBUG) << "field variable \"" << scalarFieldVariables[i]->name() << "\" has " << values[i].size() << " values: " << values[i];
  }

  // handle geometry variable
  if (geometryTable_.find(meshName) == geometryTable_.end())
  {
    std::stringstream variableName;
    variableName << "xyz";

    // communicate offset into the global values array
    int localSize = geometryFieldValues_->size()*3;
    int offset = 0;
    int globalSize = 0;

    MPI_Exscan(&localSize, &offset, 1, MPI_INT, MPI_SUM, rankSubset->mpiCommunicator());
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_INT, MPI_SUM, rankSubset->mpiCommunicator());

    LOG(DEBUG) << "define variable \"" << variableName.str() << "\", localSize: " << localSize << ", offset: " << offset << ", globalSize: " << globalSize;

    // name, global size, offset, local size
    geometryTable_[meshName] = adiosIo->DefineVariable<double>(variableName.str(), {(long unsigned int)globalSize},
                                                               {(long unsigned int)offset}, {(long unsigned int)localSize}, adios2::ConstantDims);

    if (!geometryTable_[meshName])
      LOG(ERROR) << "Failed to create geometryTable_ field variable for mesh \"" << meshName << "\".";
  }

  // handle vm variable which is the first field variable
  if (!values.empty() && vmTable_.find(meshName) == vmTable_.end())
  {
    std::stringstream variableName;
    //variableName << meshName << "_vm";
    variableName << "i";


    // communicate offset and global size
    int localSize = values[0].size();
    int offset = 0;
    int globalSize = 0;

    MPI_Exscan(&localSize, &offset, 1, MPI_INT, MPI_SUM, rankSubset->mpiCommunicator());
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_INT, MPI_SUM, rankSubset->mpiCommunicator());

    LOG(DEBUG) << "define variable \"" << variableName.str() << "\", localSize: " << localSize << ", offset: " << offset << ", globalSize: " << globalSize;

    vmTable_[meshName] = adiosIo->DefineVariable<double>(variableName.str(), {(long unsigned int)globalSize},
                                                        {(long unsigned int)offset}, {(long unsigned int)localSize}, adios2::ConstantDims);

  }

  // compute local bounding box
  if (!boundingBox.initialized)
  {
    boundingBox.min = (*geometryFieldValues_)[0];
    boundingBox.max = (*geometryFieldValues_)[0];
    boundingBox.initialized = true;
  }
  for(std::vector<Vec3>::iterator iter = geometryFieldValues_->begin(); iter != geometryFieldValues_->end(); iter++)
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

  // convert data to be send to ADIOS
  // geometryTable
  if (!geometryTableData_[meshName])
  {
    geometryTableData_[meshName] = std::make_shared<std::vector<double>>();
  }
  geometryTableData_[meshName]->resize(geometryFieldValues_->size()*3);

  // copy geometry values to buffer
  for (int i = 0; i < geometryFieldValues_->size(); i++)
  {
    for (int j = 0; j != 3; j++)
    {
      (*geometryTableData_[meshName])[3*i + j] = (*geometryFieldValues_)[i][j];
    }
  }

  // vmTable
  if (!values.empty())
  {
    // create persistent buffer, if it does not yet exist
    if (!vmTableData_[meshName])
    {
      vmTableData_[meshName] = std::make_shared<std::vector<double>>();
    }
    vmTableData_[meshName]->resize(values[0].size());

    // copy vm values to buffer
    std::copy(values[0].begin(), values[0].end(), vmTableData_[meshName]->begin());

    assert(vmTable_[meshName]);
    adiosWriter->Put<double>(vmTable_[meshName], vmTableData_[meshName]->data());  // adios2::Mode::Sync
  }

  assert(geometryTable_[meshName]);
  adiosWriter->Put<double>(geometryTable_[meshName], geometryTableData_[meshName]->data());

}
#endif

}  // namespace
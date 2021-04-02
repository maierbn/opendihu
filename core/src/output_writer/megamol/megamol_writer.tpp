#include "output_writer/megamol/megamol_writer.h"

#include "output_writer/megamol/loop_collect_field_variables.h"

#ifdef HAVE_ADIOS
#include <adios2.h>
#endif

namespace OutputWriter
{

#ifdef HAVE_ADIOS
template<typename FunctionSpaceType, typename FieldVariablesForOutputWriterType>
void MegaMolWriter<FunctionSpaceType, FieldVariablesForOutputWriterType>::
outputData(FieldVariablesForOutputWriterType fieldVariables, std::string meshName, std::shared_ptr<FunctionSpaceType> functionSpace,
           PythonConfig specificSettings, MegaMolWriterContext &megaMolWriterContext)
{
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField;
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> scalarFieldVariables;

  // collect the geometryField and all scalar field variables for the current mesh
  MegaMolLoopOverTuple::loopCollectFieldVariables(fieldVariables, meshName, geometryField, scalarFieldVariables);

  // retrieve the new geometry field values
  int oldSize = megaMolWriterContext.geometryFieldVectorValues.size();
  geometryField->getValuesWithoutGhosts(megaMolWriterContext.geometryFieldVectorValues);

  LOG(DEBUG) << "nPointsPerCoordinateDirection: " << megaMolWriterContext.nPointsPerCoordinateDirection;
  LOG(DEBUG) << "number of scalar field variables: " << scalarFieldVariables.size();

  // get first other scalar field variable
  if (scalarFieldVariables.size() > 0 && scalarFieldVariables.size() <= 1)
  {
    scalarFieldVariables[0]->getValuesWithoutGhosts(megaMolWriterContext.vmValues);
    LOG(DEBUG) << "mesh \"" << meshName << "\", retrieve field variable \"" << scalarFieldVariables[0]->name() << "\", n values: " << megaMolWriterContext.vmValues.size();

    megaMolWriterContext.nPointsPerCoordinateDirection[0] = geometryField->functionSpace()->meshPartition()->nNodesGlobal(0);
    megaMolWriterContext.nPointsPerCoordinateDirection[1] = 0;
    megaMolWriterContext.nPointsPerCoordinateDirection[2] = 0;
  }
  else
  {
    megaMolWriterContext.nPointsPerCoordinateDirection[0] = geometryField->functionSpace()->meshPartition()->nNodesGlobal(0);

    if (geometryField->functionSpace()->dim() >= 2)
      megaMolWriterContext.nPointsPerCoordinateDirection[1] = geometryField->functionSpace()->meshPartition()->nNodesGlobal(1);
    else
      megaMolWriterContext.nPointsPerCoordinateDirection[1] = 1;

    if (geometryField->functionSpace()->dim() >= 3)
      megaMolWriterContext.nPointsPerCoordinateDirection[2] = geometryField->functionSpace()->meshPartition()->nNodesGlobal(2);
    else
      megaMolWriterContext.nPointsPerCoordinateDirection[1] = 1;

    LOG(DEBUG) << "evaluate names of field variables";
    for (typename std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>>::iterator iter = scalarFieldVariables.begin();
         iter != scalarFieldVariables.end(); iter++)
    {
      LOG(DEBUG) << "field variable \"" << (*iter)->name() << "\".";
      if ((*iter)->name() == "Vm")
      {
        (*iter)->getValuesWithoutGhosts(megaMolWriterContext.vmValues);
        LOG(DEBUG) << "->vmValues";
        LOG(DEBUG) << "n values: " << megaMolWriterContext.vmValues.size();
      }
      else if ((*iter)->name() == "phi_e")
      {
        (*iter)->getValuesWithoutGhosts(megaMolWriterContext.emgValues);
        LOG(DEBUG) << "->phi_e";
        LOG(DEBUG) << "n values: " << megaMolWriterContext.emgValues.size();
      }
      else if ((*iter)->name() == "transmembraneFlow")
      {
        (*iter)->getValuesWithoutGhosts(megaMolWriterContext.transmembraneFlowValues);
        LOG(DEBUG) << "->trans";
        LOG(DEBUG) << "n values: " << megaMolWriterContext.transmembraneFlowValues.size();
      }
    }
    //LOG(INFO) << "mesh \"" << meshName << "\", retrieve field variable \"" << scalarFieldVariables[0]->name() << "\", n values: " << megaMolWriterContext.scalarFieldVariableValues.size();
  }

  int nNodes = megaMolWriterContext.geometryFieldVectorValues.size();
  megaMolWriterContext.rankNo.resize(nNodes, DihuContext::ownRankNoCommWorld());

  LOG(DEBUG) << "MegaMolWriter::outputData, mesh \"" << meshName << "\".";

  if (oldSize == 0)
    return;

  // compute distance between the last points of the previous and the new fiber
  double approximateDistanceBetweenFibers = MathUtility::distance<3>(
    megaMolWriterContext.geometryFieldVectorValues[oldSize-1],
    megaMolWriterContext.geometryFieldVectorValues[megaMolWriterContext.geometryFieldVectorValues.size()-1]
  );

  // store value as minimum to approximateDistanceBetweenFibers
  if (megaMolWriterContext.approximateDistanceBetweenFibers == -1)
  {
    megaMolWriterContext.approximateDistanceBetweenFibers = approximateDistanceBetweenFibers;
  }
  else
  {
    megaMolWriterContext.approximateDistanceBetweenFibers = std::min(megaMolWriterContext.approximateDistanceBetweenFibers, approximateDistanceBetweenFibers);
  }
}
#endif

}  // namespace

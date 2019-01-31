#include "output_writer/megamol/megamol_writer.h"

#include "output_writer/megamol/loop_collect_field_variables.h"

#ifdef HAVE_ADIOS
#include <adios2.h>
#endif

namespace OutputWriter
{

#ifdef HAVE_ADIOS
template<typename FunctionSpaceType, typename OutputFieldVariablesType>
void MegaMolWriter<FunctionSpaceType, OutputFieldVariablesType>::
outputData(OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<FunctionSpaceType> functionSpace,
           PythonConfig specificSettings, MegaMolWriterContext &megaMolWriterContext)
{
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField;
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> scalarFieldVariables;

  // collect the geometryField and all scalar field variables for the current mesh
  MegaMolLoopOverTuple::loopCollectFieldVariables(fieldVariables, meshName, geometryField, scalarFieldVariables);

  // retrieve the new geometry field values
  int oldSize = megaMolWriterContext.geometryFieldValues.size();
  geometryField->getValuesWithoutGhosts(megaMolWriterContext.geometryFieldValues);

  // get first other scalar field variable
  if (scalarFieldVariables.size() > 0)
  {
    LOG(DEBUG) << "mesh \"" << meshName << "\", retrieve field variable \"" << scalarFieldVariables[0]->name() << "\"";
    scalarFieldVariables[0]->getValuesWithoutGhosts(megaMolWriterContext.scalarFieldVariableValues);
  }

  LOG(DEBUG) << "MegaMolWriter::outputData, mesh \"" << meshName << "\".";

  if (oldSize == 0)
    return;

  // compute distance between the last points of the previous and the new fiber
  double approximateDistanceBetweenFibers = MathUtility::distance<3>(
    megaMolWriterContext.geometryFieldValues[oldSize-1],
    megaMolWriterContext.geometryFieldValues[megaMolWriterContext.geometryFieldValues.size()-1]
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

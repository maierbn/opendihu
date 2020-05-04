#include "mesh/mapping_between_meshes/manager/01_manager_implementation.h"

#include <memory>

#include "easylogging++.h"
#include "mesh/structured_regular_fixed.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "utility/vector_operators.h"
#include "partition/partitioned_petsc_vec/values_representation.h"

namespace MappingBetweenMeshes
{

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
std::shared_ptr<MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>> ManagerImplementation::
createMappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget)
{
  std::string sourceMeshName = functionSpaceSource->meshName();
  std::string targetMeshName = functionSpaceTarget->meshName();

  // check if the mapping already existed
  if (this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].mapping)
  {
    LOG(WARNING) << "Mapping from mesh \"" << sourceMeshName << "\" to mesh \"" << targetMeshName << "\" is already defined.";
  }

  // get options for the mapping
  double xiTolerance = this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].xiTolerance;
  bool enableWarnings = this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].enableWarnings;
  bool compositeUseOnlyInitializedMappings = this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].compositeUseOnlyInitializedMappings;
  bool isEnabledFixUnmappedDofs = this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].isEnabledFixUnmappedDofs;

  LOG(INFO) << "create MappingBetweenMeshes \"" << sourceMeshName << "\" (" << FunctionSpaceSourceType::dim() << "D, " << functionSpaceSource->nNodesGlobal() << " nodes) -> \""
     << targetMeshName << "\" (" << FunctionSpaceTargetType::dim() << "D, " << functionSpaceTarget->nNodesGlobal() << " nodes), xiTolerance: " << xiTolerance;

  // create the mapping under the given source and target mesh names
  this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].mapping = std::static_pointer_cast<MappingBetweenMeshesBase>(
    std::make_shared<MappingBetweenMeshes<FunctionSpaceSourceType,FunctionSpaceTargetType>>(functionSpaceSource, functionSpaceTarget,
                                                                                            xiTolerance, enableWarnings, compositeUseOnlyInitializedMappings, isEnabledFixUnmappedDofs)
  );

  // log event, to be included in the log file
  addLogEntryMapping(functionSpaceSource, functionSpaceTarget, mappingLogEntry_t::logEvent_t::eventCreateMapping);

  return std::static_pointer_cast<MappingBetweenMeshes<FunctionSpaceSourceType,FunctionSpaceTargetType>>(
    this->mappingsBetweenMeshes_[sourceMeshName][targetMeshName].mapping);
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
std::shared_ptr<MappingBetweenMeshes<typename FunctionSpaceSourceType::FunctionSpace, typename FunctionSpaceTargetType::FunctionSpace>> ManagerImplementation::
mappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                     std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget)
{
  typedef MappingBetweenMeshes<typename FunctionSpaceSourceType::FunctionSpace, typename FunctionSpaceTargetType::FunctionSpace> MappingType;

  assert(functionSpaceSource);
  assert(functionSpaceTarget);

  // get mesh names
  std::string sourceMeshName = functionSpaceSource->meshName();
  std::string targetMeshName = functionSpaceTarget->meshName();

  // check if the mapping exists already
  if (mappingsBetweenMeshes_.find(sourceMeshName) != mappingsBetweenMeshes_.end())
  {
    if (mappingsBetweenMeshes_[sourceMeshName].find(targetMeshName) != mappingsBetweenMeshes_[sourceMeshName].end())
    {
      std::shared_ptr<MappingBetweenMeshesBase> mappingBase = mappingsBetweenMeshes_[sourceMeshName][targetMeshName].mapping;

      // if mapping has already been created, return it
      if (mappingBase)
      {
        return std::static_pointer_cast<MappingType>(mappingBase);
      }
    }
  }

  // mapping between meshes has not yet been created, create now
  // if it does not yet exist, output message and create it
  LOG(DEBUG) << "Mapping from mesh \"" << sourceMeshName << "\" to \"" << targetMeshName
    << "\" was not initialized. Initializing now. Specify MappingsBetweenMeshes { \"" << sourceMeshName << "\" : \"" << targetMeshName << "\" } as top level object of the python config. "
    << "(It could be that this was done, but because of MultipleInstances in the OperatorSplitting or Coupling, only the first mesh mapping got initialized.)";

  // create the mapping
  std::shared_ptr<MappingType> mapping;
  mapping = createMappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>(
    functionSpaceSource, functionSpaceTarget
  );
  return mapping;
}

//! prepare a mapping to the fieldVariableTarget, this zeros the targetFactorSum for the field variable
template<typename FieldVariableTargetType>
void ManagerImplementation::
prepareMappingLowToHigh(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget)
{
  Control::PerformanceMeasurement::start("durationMapPrepare");

  VLOG(1) << "prepareMappingLowToHigh, fieldVariableTarget: " << fieldVariableTarget->name() << " componentNoTarget: " << componentNoTarget;

  std::string fieldVariableTargetName = fieldVariableTarget->name();
  // std::map<std::string, std::shared_ptr<FieldVariable::FieldVariableBase>> targetFactorSum_;

  // if the targetFactorSum field variable does not yet exist, create it
  std::string targetFactorSumName = fieldVariableTarget->functionSpace()->meshName()+std::string("_")+fieldVariableTarget->name();
  if (targetFactorSum_.find(targetFactorSumName) == targetFactorSum_.end())
  {
    std::stringstream name;
    name << "targetFactorSum_" << targetFactorSumName;

    std::vector<std::string> componentNames(1, "0");

    targetFactorSum_[targetFactorSumName] = std::static_pointer_cast<FieldVariable::FieldVariableBase>(
      std::make_shared<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>>(*fieldVariableTarget, name.str(), componentNames)
    );
  }

  std::shared_ptr<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>> targetFactorSum
    = std::static_pointer_cast<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>>(targetFactorSum_[targetFactorSumName]);

  VLOG(1) << "prepareMappingLowToHigh, set all entries of " << targetFactorSum->name() << " and " << fieldVariableTarget->name() << " to 0";

  // set all entries to 0
  targetFactorSum->zeroEntries();
  targetFactorSum->zeroGhostBuffer();

  // zero the entries of the component that will be set
  zeroTargetFieldVariable(fieldVariableTarget, componentNoTarget);

  Control::PerformanceMeasurement::stop("durationMapPrepare");
}

template<typename FieldVariableTargetType>
void ManagerImplementation::
zeroTargetFieldVariable(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget)
{
  // zero the entries of the component that will be set
  if (componentNoTarget == -1)
  {
    fieldVariableTarget->zeroEntries();
  }
  else 
  {
    int nDofsLocalWithoutGhosts = fieldVariableTarget->functionSpace()->nDofsLocalWithoutGhosts();
    std::vector<double> zeros(nDofsLocalWithoutGhosts, 0);
    fieldVariableTarget->setValuesWithoutGhosts(componentNoTarget, zeros, INSERT_VALUES);
  }
  fieldVariableTarget->zeroGhostBuffer();
}

// helper function, calls the map function of the mapping if field variables have same number of components
template<typename FieldVariableSourceType, typename FieldVariableTargetType, typename Dummy = FieldVariableSourceType>
struct MapLowToHighDimensionAllComponents
{

  // helper function, does nothing
  template<typename T1, typename T2>
  static void call(
    T1 mapping,
    std::shared_ptr<FieldVariableSourceType> fieldVariableSource, std::shared_ptr<FieldVariableTargetType> fieldVariableTarget,
    T2 targetFactorSum)
  {
    LOG(FATAL) << "Number of components of field variables does not match"
      << "(" << FieldVariableSourceType::nComponents() << " != " << FieldVariableTargetType::nComponents() << "),"
      << " but a mapping between all components was requested.";
  }
};

// helper function, calls the map function of the mapping if field variables have same number of components
template<typename FieldVariableSourceType, typename FieldVariableTargetType>
struct MapLowToHighDimensionAllComponents<FieldVariableSourceType, FieldVariableTargetType,
    typename std::enable_if<FieldVariableSourceType::nComponents() == FieldVariableTargetType::nComponents(),FieldVariableSourceType>::type>
{

  // actual function
  static void call(
    std::shared_ptr<MappingBetweenMeshes<typename FieldVariableSourceType::FunctionSpace, typename FieldVariableTargetType::FunctionSpace>> mapping,
    std::shared_ptr<FieldVariableSourceType> fieldVariableSource, std::shared_ptr<FieldVariableTargetType> fieldVariableTarget,
    std::shared_ptr<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>> targetFactorSum)
  {
    mapping->template mapLowToHighDimension<FieldVariableSourceType::nComponents()>(
      *fieldVariableSource, *fieldVariableTarget, *targetFactorSum
    );
  }

};

//! map data from the source to the target field variable. This has to be called between prepareMapping and finalizeMapping, can be called multiple times with different source meshes.
template<typename FieldVariableSourceType, typename FieldVariableTargetType>
void ManagerImplementation::
mapLowToHighDimension(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
                      std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget)
{
  // source = lower dimension
  // target = higher dimension

  std::string sourceMeshName = fieldVariableSource->functionSpace()->meshName();
  std::string targetMeshName = fieldVariableTarget->functionSpace()->meshName();

  typedef MappingBetweenMeshes<typename FieldVariableSourceType::FunctionSpace, typename FieldVariableTargetType::FunctionSpace> MappingType;

  std::shared_ptr<MappingType> mapping = this->mappingBetweenMeshes<typename FieldVariableSourceType::FunctionSpace, typename FieldVariableTargetType::FunctionSpace>(
    fieldVariableSource->functionSpace(), fieldVariableTarget->functionSpace()
  );

  Control::PerformanceMeasurement::start("durationMap");

  // assert that targetFactorSum_ field variable exists, this should have been created by prepareMapping()
  std::string targetFactorSumName = targetMeshName+std::string("_")+fieldVariableTarget->name();
  assert(targetFactorSum_.find(targetFactorSumName) != targetFactorSum_.end());

  std::shared_ptr<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>> targetFactorSum
    = std::static_pointer_cast<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>>(
      targetFactorSum_[targetFactorSumName]);

  // assert that both or none of the componentNos are -1
  assert((componentNoSource == -1) == (componentNoTarget == -1));

  // if all components should be transferred
  if (componentNoSource == -1 && componentNoTarget == -1)
  {
    LOG(DEBUG) << "map low to high dimension, " << sourceMeshName << " -> " << targetMeshName << ", all components ";
    LOG(DEBUG) << "mapping: " << mapping;
    LOG(DEBUG) << "fieldVariableSource: " << fieldVariableSource;
    LOG(DEBUG) << "fieldVariableTarget: " << fieldVariableTarget;

    // call the method of the mapping that does the actual data transfer
    MapLowToHighDimensionAllComponents<FieldVariableSourceType,FieldVariableTargetType>::call(mapping, fieldVariableSource, fieldVariableTarget, targetFactorSum);
    //mapping->template mapLowToHighDimension<FieldVariableSourceType::nComponents()>(
    //  *fieldVariableSource, *fieldVariableTarget, *targetFactorSum
    //);
  }
  else
  {
    // if only the specified components should be transferred

    // call the method of the mapping that does the actual data transfer
    mapping->template mapLowToHighDimension<FieldVariableSourceType::nComponents(),FieldVariableTargetType::nComponents()>(
      *fieldVariableSource, componentNoSource, *fieldVariableTarget, componentNoTarget, *targetFactorSum
    );
  }

  Control::PerformanceMeasurement::stop("durationMap");
}

// helper function, calls the map function of the mapping if field variables have same number of components
template<typename FieldVariableSourceType, typename FieldVariableTargetType, typename Dummy=void>
struct MapHighToLowDimensionAllComponents
{
  // helper function, does nothing
  template<typename T1>
  static void call(
    T1 mapping,
    std::shared_ptr<FieldVariableSourceType> fieldVariableSource, std::shared_ptr<FieldVariableTargetType> fieldVariableTarget)
  {
    LOG(FATAL) << "Number of components of field variables does not match"
      << "(" << FieldVariableSourceType::nComponents() << " != " << FieldVariableTargetType::nComponents() << "),"
      << " but a mapping between all components was requested.";
  }
};

// helper function, calls the map function of the mapping if field variables have same number of components
template<typename FieldVariableSourceType, typename FieldVariableTargetType>
struct MapHighToLowDimensionAllComponents<FieldVariableSourceType,FieldVariableTargetType,
  typename std::enable_if<FieldVariableSourceType::nComponents() == FieldVariableTargetType::nComponents(),void>::type
>
{
  static void call(
    std::shared_ptr<MappingBetweenMeshes<typename FieldVariableTargetType::FunctionSpace, typename FieldVariableSourceType::FunctionSpace>> mapping,
    std::shared_ptr<FieldVariableSourceType> fieldVariableSource, std::shared_ptr<FieldVariableTargetType> fieldVariableTarget)
  {
    mapping->template mapHighToLowDimension<FieldVariableSourceType::nComponents()>(
      *fieldVariableSource, *fieldVariableTarget
    );
  }
};

//! map specific componentNo. This has to be called between prepareMapping and finalizeMapping, can be called multiple times with different source meshes.
template<typename FieldVariableSourceType, typename FieldVariableTargetType>
void ManagerImplementation::
mapHighToLowDimension(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
                      std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget)
{
  // source = higher dimension
  // target = lower dimension

  std::string sourceMeshName = fieldVariableSource->functionSpace()->meshName();
  std::string targetMeshName = fieldVariableTarget->functionSpace()->meshName();

  // here the mapping is defined with first FunctionSpace being the target function space and second FunctionSpace being the source function space.
  typedef MappingBetweenMeshes<typename FieldVariableTargetType::FunctionSpace, typename FieldVariableSourceType::FunctionSpace> MappingType;

  std::shared_ptr<MappingType> mapping = this->mappingBetweenMeshes<typename FieldVariableTargetType::FunctionSpace, typename FieldVariableSourceType::FunctionSpace>(
    fieldVariableTarget->functionSpace(), fieldVariableSource->functionSpace()
  );

  Control::PerformanceMeasurement::start("durationMap");

  // assert that both or none of the componentNos are -1
  assert((componentNoSource == -1) == (componentNoTarget == -1));

  // call the method of the mapping that does the actual data transfer
  if (componentNoSource == -1 && componentNoTarget == -1)
  {
    MapHighToLowDimensionAllComponents<FieldVariableSourceType,FieldVariableTargetType>::call(mapping, fieldVariableSource, fieldVariableTarget);
  }
  else
  {
    mapping->template mapHighToLowDimension<FieldVariableSourceType::nComponents()>(
      *fieldVariableSource, componentNoSource, *fieldVariableTarget, componentNoTarget
    );
  }

  Control::PerformanceMeasurement::stop("durationMap");
}

//! finalize the mapping to the fieldVariableTarget, this computes the final values at the dofs from the accumulated values by dividing by the targetFactorSums
template<typename FieldVariableTargetType>
void ManagerImplementation::
finalizeMappingLowToHigh(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget)
{
  if (componentNoTarget == -1)
  {
    finalizeMappingLowToHigh(fieldVariableTarget);
    return;
  }

  VLOG(1) << "finalizeMappingLowToHigh, fieldVariableTarget: " << fieldVariableTarget->name() << ", componentNoTarget: " << componentNoTarget;

  Control::PerformanceMeasurement::start("durationMapFinalize");

  std::string targetFactorSumName = fieldVariableTarget->functionSpace()->meshName()+std::string("_")+fieldVariableTarget->name();
  // assert that targetFactorSum_ field variable exists, this should have been created by prepareMapping()
  assert(targetFactorSum_.find(targetFactorSumName) != targetFactorSum_.end());


  std::shared_ptr<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>> targetFactorSum
    = std::static_pointer_cast<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>>(targetFactorSum_[targetFactorSumName]);

  const dof_no_t nDofsLocalTarget = fieldVariableTarget->nDofsLocalWithoutGhosts();

  // communicate ghost values
  fieldVariableTarget->finishGhostManipulation();
  targetFactorSum->finishGhostManipulation();

  // compute final values by dividing by factorSum at each target dof
  std::vector<double> targetValues;
  fieldVariableTarget->getValuesWithoutGhosts(componentNoTarget, targetValues);

  std::vector<double> targetFactorSums;
  targetFactorSum->getValuesWithoutGhosts(targetFactorSums);

  std::stringstream info;

  for (dof_no_t targetDofNoLocal = 0; targetDofNoLocal != nDofsLocalTarget; targetDofNoLocal++)
  {
    info << "  target dof " << targetDofNoLocal << ", divide value " << targetValues[targetDofNoLocal]
      << " by " << targetFactorSums[targetDofNoLocal] << ": " << targetValues[targetDofNoLocal]/targetFactorSums[targetDofNoLocal];


    VLOG(2) << "  target dof " << targetDofNoLocal << ", divide value " << targetValues[targetDofNoLocal]
      << " by " << targetFactorSums[targetDofNoLocal] << ": " << targetValues[targetDofNoLocal]/targetFactorSums[targetDofNoLocal];
    if (fabs(targetFactorSums[targetDofNoLocal]) > 1e-12)
    {
      targetValues[targetDofNoLocal] /= targetFactorSums[targetDofNoLocal];
    }
    else
    {
#ifndef NDEBUG
      // output warning message, compute helper variables

      // get node no from dof no, this is needed for the global coordinates later
      //node_no_t nodeNoLocal = int(targetDofNoLocal/fieldVariableTarget->functionSpace()->nDofsPerNode());

      // get current global node coordinates
      std::stringstream s;
      /*for (int i = 0; i < FieldVariableTargetType::FunctionSpace::dim(); i++)
      {
        if (i != 0)
        {
          s << ",";
        }
        s << fieldVariableTarget->functionSpace()->meshPartition()->nNodesGlobal(i);
      }*/

      // output the warning
      LOG(WARNING) << "In mapping to " << fieldVariableTarget->name() << "." << componentNoTarget
        << " (" << fieldVariableTarget->functionSpace()->meshName() << "), no values for target dof " << targetDofNoLocal
        //<< ", coordinates global: " << fieldVariableTarget->functionSpace()->meshPartition()->getCoordinatesGlobal(nodeNoLocal)
        //<< " of (" << s.str() << ")! "
        << ". Assuming 0.0.";
#endif
    }
  }

  VLOG(1) << "low to high, set targetValues: " << targetValues;

  // set the computed values
  fieldVariableTarget->setValuesWithoutGhosts(componentNoTarget, targetValues);

  Control::PerformanceMeasurement::stop("durationMapFinalize");
}

//! finalize the mapping to the fieldVariableTarget, this computes the final values at the dofs from the accumulated values by dividing by the targetFactorSums
template<typename FieldVariableTargetType>
void ManagerImplementation::
finalizeMappingLowToHigh(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget)
{
  Control::PerformanceMeasurement::start("durationMapFinalize");

  VLOG(1) << "finalizeMappingLowToHigh, fieldVariableTarget: " << fieldVariableTarget->name();

  // assert that targetFactorSum_ field variable exists, this should have been created by prepareMapping()
  std::string targetFactorSumName = fieldVariableTarget->functionSpace()->meshName()+std::string("_")+fieldVariableTarget->name();
  assert(targetFactorSum_.find(targetFactorSumName) != targetFactorSum_.end());


  std::shared_ptr<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>> targetFactorSum
    = std::static_pointer_cast<FieldVariable::FieldVariable<typename FieldVariableTargetType::FunctionSpace,1>>(targetFactorSum_[targetFactorSumName]);

  const dof_no_t nDofsLocalTarget = fieldVariableTarget->nDofsLocalWithoutGhosts();

  // communicate ghost values
  fieldVariableTarget->finishGhostManipulation();
  targetFactorSum->finishGhostManipulation();

  // compute final values by dividing by factorSum
  std::vector<VecD<FieldVariableTargetType::nComponents()>> targetValues;
  fieldVariableTarget->getValuesWithoutGhosts(targetValues);

  std::vector<double> targetFactorSums;
  targetFactorSum->getValuesWithoutGhosts(targetFactorSums);

  for (dof_no_t targetDofNoLocal = 0; targetDofNoLocal != nDofsLocalTarget; targetDofNoLocal++)
  {
    VLOG(2) << "  target dof " << targetDofNoLocal << ", divide value " << targetValues[targetDofNoLocal]
      << " by " << targetFactorSums[targetDofNoLocal] << ": " << targetValues[targetDofNoLocal]/targetFactorSums[targetDofNoLocal];
    if (fabs(targetFactorSums[targetDofNoLocal]) > 1e-12)
    {
      targetValues[targetDofNoLocal] /= targetFactorSums[targetDofNoLocal];
    }
    else
    {

#ifndef NDEBUG
      // output warning message, compute helper variables

      // get node no from dof no, this is needed for the global coordinates later
      //node_no_t nodeNoLocal = int(targetDofNoLocal/fieldVariableTarget->functionSpace()->nDofsPerNode());

      // get current global node coordinates
      std::stringstream s;
/*      for (int i = 0; i < FieldVariableTargetType::FunctionSpace::dim(); i++)
      {
        if (i != 0)
        {
          s << ",";
        }
        s << fieldVariableTarget->functionSpace()->meshPartition()->nNodesGlobal(i);
      }
*/
      // output the warning
      LOG(WARNING) << "In mapping to " << fieldVariableTarget->name()
        << " (" << fieldVariableTarget->functionSpace()->meshName() << "), no values for target dof " << targetDofNoLocal
        //<< ", coordinates global: " << fieldVariableTarget->functionSpace()->meshPartition()->getCoordinatesGlobal(nodeNoLocal)
        //<< " of (" << s.str() << ")! Assuming 0.0.";
        << ". Assuming 0.0.";
#endif
    }
  }

  VLOG(1) << "low to high all, set targetValues: " << targetValues;

  // set the computed values
  fieldVariableTarget->setValuesWithoutGhosts(targetValues);

  Control::PerformanceMeasurement::stop("durationMapFinalize");
}

}   // namespace

#include "mesh/mapping_between_meshes/manager/03_manager_implementation.h"

#include <memory>

#include "easylogging++.h"
#include "mesh/structured_regular_fixed.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "utility/vector_operators.h"
#include "partition/partitioned_petsc_vec/values_representation.h"

namespace MappingBetweenMeshes
{

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

  std::string targetMeshName = fieldVariableTarget->functionSpace()->meshName();
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
      // if there was a default value specified, set the target value to the default value
      if (defaultValues_.find(targetMeshName) != defaultValues_.end())
      {
        // set to default value if there is one
        targetValues[targetDofNoLocal] = defaultValues_[targetMeshName];
      }

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

  std::string targetMeshName = fieldVariableTarget->functionSpace()->meshName();

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
      // if there was a default value specified, set the target value to the default value
      if (defaultValues_.find(targetMeshName) != defaultValues_.end())
      {
        // set to default value if there is one
        targetValues[targetDofNoLocal].fill(defaultValues_[targetMeshName]);
      }

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

//! repair invalid geometry values for 1D fibers after mapping from a 3D mesh
template<typename FieldVariableTargetType>
void ManagerImplementation::
repairMappedGeometryFibers(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget)
{
  // this method fixes the geometry of fibers if single points have wrong values after the mapping
  // this occurs only for highly irregularly shaped elements
  LOG(DEBUG) << "repairMappedGeometryFibers " << fieldVariableTarget->name() << " on " << fieldVariableTarget->functionSpace()->meshName();
  
  const dof_no_t nDofsLocalTarget = fieldVariableTarget->nDofsLocalWithoutGhosts();

  using Vec = VecD<FieldVariableTargetType::nComponents()>;
  const int D = FieldVariableTargetType::nComponents();

  std::vector<Vec> targetValues;
  fieldVariableTarget->getValuesWithoutGhosts(targetValues);

  // compute median distance between points and compute median position
  static std::set<double> distances;
  static std::set<double> xpos;
  static std::set<double> ypos;
  static std::set<double> zpos;
  
  Vec3 sum = MathUtility::transformToD<3,D>(targetValues[0]);
  // iterate over all elements of the fiber
  for (dof_no_t targetDofNoLocal = 1; targetDofNoLocal != nDofsLocalTarget; targetDofNoLocal++)
  {
    Vec3 position = MathUtility::transformToD<3,D>(targetValues[targetDofNoLocal]);
    Vec3 previousPosition = MathUtility::transformToD<3,D>(targetValues[targetDofNoLocal-1]);
    
    double distance = MathUtility::distance<3>(position,previousPosition);
    distances.insert(distance);
    xpos.insert(position[0]);
    ypos.insert(position[1]);
    zpos.insert(position[2]);
  }
  
  // determine the median element length (distance)
  std::set<double>::const_iterator distanceIter = distances.begin();
  int nEntries = distances.size();
  for (int i = 0; i < nEntries/2; i++, distanceIter++);
  double medianDistance = *distanceIter;
  
  // determine the median x position
  std::set<double>::const_iterator xposIter = xpos.begin();
  for (int i = 0; i < xpos.size()/2; i++, xposIter++);
  double medianXPos = *xposIter;
  
  // determine the median y position
  std::set<double>::const_iterator yposIter = ypos.begin();
  for (int i = 0; i < ypos.size()/2; i++, yposIter++);
  double medianYPos = *yposIter;
  
  // determine the median z position
  std::set<double>::const_iterator zposIter = zpos.begin();
  for (int i = 0; i < zpos.size()/2; i++, zposIter++);
  double medianZPos = *zposIter;
  
  LOG(DEBUG) << "distances: " << distances;
  LOG(DEBUG) << "medianDistance: " << medianDistance;
  
  // check if any distances are invalid and fix them
  
  dof_no_t lastValidDof = -1;
  dof_no_t validDofFiberBegin = -1;
  
  // loop over all nodes
  for (dof_no_t targetDofNoLocal = 0; targetDofNoLocal != nDofsLocalTarget; targetDofNoLocal++)
  {
    Vec3 position = MathUtility::transformToD<3,D>(targetValues[targetDofNoLocal]);
    
    // if point is invalid (position is for away from median in all dimensions)
    if (fabs(position[0] - medianXPos) > 5 && fabs(position[1] - medianYPos) > 5 && fabs(position[2] - medianZPos) > 5)
    {
      // do nothing, position will get fixed at next valid dof
    }
    else
    {
      // if point is valid
      
      // if there is at least one valid point beforehand
      if (lastValidDof != -1)
      {
        // fix all previous invalid points
        Vec3 lastValidPosition = MathUtility::transformToD<3,D>(targetValues[lastValidDof]);
        Vec3 v = -lastValidPosition + position;
        
        // loop over all previous invalid points
        int index = 1;
        int nInvalidElements = targetDofNoLocal - lastValidDof;
        for (dof_no_t dofNoLocal = lastValidDof+1; dofNoLocal < targetDofNoLocal; dofNoLocal++, index++)
        {
          Vec3 fixedPosition = lastValidPosition + v * index / nInvalidElements;
          
          LOG(DEBUG) << "Fix dof " << dofNoLocal << "/" << nDofsLocalTarget << ": " 
            << fixedPosition << ", previous: " << targetValues[dofNoLocal];
            
          // assign fixed position to target values entry
          for (int i = 0; i < 3; i++)
            targetValues[dofNoLocal][i] = fixedPosition[i];
        }
      }
        
      // store valid dof no
      lastValidDof = targetDofNoLocal;
      if (validDofFiberBegin == -1)
        validDofFiberBegin = targetDofNoLocal;
    }
  }
  
  // fix invalid points at the beginning of the fiber
  if (validDofFiberBegin > 0 && validDofFiberBegin + 1 < nDofsLocalTarget)
  {
    // beginning of fiber (i=invalid, validDofFiberBegin=4):
    // i i i i p0 p1 p2 ...
    
    Vec3 p0 = MathUtility::transformToD<3,D>(targetValues[validDofFiberBegin]);
    Vec3 p1 = MathUtility::transformToD<3,D>(targetValues[validDofFiberBegin+1]);
    Vec3 v = -p1 + p0;
    
    // loop over the first invalid points until the first valid
    for (dof_no_t dofNoLocal = 0; dofNoLocal < validDofFiberBegin; dofNoLocal++)
    {
      Vec3 fixedPosition = p0 + v * (validDofFiberBegin - dofNoLocal);
      
      LOG(DEBUG) << "Fix dof " << dofNoLocal << "/" << nDofsLocalTarget << " (at beginning): " 
        << fixedPosition << ", previous: " << targetValues[dofNoLocal];
        
      // assign fixed position to target values entry
      for (int i = 0; i < 3; i++)
        targetValues[dofNoLocal][i] = fixedPosition[i];
    }
  }
  
  // fix invalid points at the end of the fiber
  if (lastValidDof < nDofsLocalTarget-1 && lastValidDof > 0)
  {
    // end of fiber (i=invalid, lastValidDof=nDofsLocalTarget-5):
    // ... p0 p1 i i i i
    
    Vec3 p0 = MathUtility::transformToD<3,D>(targetValues[lastValidDof-1]);
    Vec3 p1 = MathUtility::transformToD<3,D>(targetValues[lastValidDof]);
    Vec3 v = -p0 + p1;
    
    // loop over the invalid points at the end
    for (dof_no_t dofNoLocal = lastValidDof+1; dofNoLocal < nDofsLocalTarget; dofNoLocal++)
    {
      Vec3 fixedPosition = p1 + v * (dofNoLocal - lastValidDof);
      
      LOG(DEBUG) << "Fix dof " << dofNoLocal << "/" << nDofsLocalTarget << " (at end): " 
        << fixedPosition << ", previous: " << targetValues[dofNoLocal];
      
      // assign fixed position to target values entry
      for (int i = 0; i < 3; i++)
        targetValues[dofNoLocal][i] = fixedPosition[i];
    }
  }
  
  // set the computed values
  fieldVariableTarget->setValuesWithoutGhosts(targetValues);
}

}   // namespace

#include "mesh/mapping_between_meshes/mapping/00_construct.h"

#include "control/diagnostic_tool/performance_measurement.h"

#include "utility/vector_operators.h"
#include "control/dihu_context.h"
#include "mesh/type_traits.h"
#include "mesh/mapping_between_meshes/manager/04_manager.h"
#include "mesh/mapping_between_meshes/manager/target_element_no_estimator.h"

namespace MappingBetweenMeshes
{

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
MappingBetweenMeshesConstruct<FunctionSpaceSourceType, FunctionSpaceTargetType>::
MappingBetweenMeshesConstruct(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                              std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget,
                              double xiTolerance, bool enableWarnings, bool compositeUseOnlyInitializedMappings,
                              bool isEnabledFixUnmappedDofs) :
  functionSpaceSource_(functionSpaceSource),
  functionSpaceTarget_(functionSpaceTarget)
{
  // for composite meshes if the option compositeUseOnlyInitializedMappings is set, do not create the mapping here
  if (Mesh::isComposite<std::shared_ptr<FunctionSpaceSourceType>>::value && compositeUseOnlyInitializedMappings)
  {
    LOG(DEBUG) << "Do not initialize mapping here, it will be done in the child class from existing mappings of the sub meshes of the composite mesh.";
  }
  else 
  {
    // create the mapping
    Control::PerformanceMeasurement::start("durationComputeMappingBetweenMeshes");

    const dof_no_t nDofsLocalSource = functionSpaceSource->nDofsLocalWithoutGhosts();
    const dof_no_t nDofsLocalTarget = functionSpaceTarget->nDofsLocalWithoutGhosts();
    const int nDofsPerTargetElement = FunctionSpaceTargetType::nDofsPerElement();

    LOG(DEBUG) << "create MappingBetweenMeshes \"" << functionSpaceSource->meshName() << "\" and \""
      << functionSpaceTarget->meshName() << "\" xiTolerance: " << xiTolerance << ", is composite: " << Mesh::isComposite<std::shared_ptr<FunctionSpaceSourceType>>::value
      << ", compositeUseOnlyInitializedMappings: " << compositeUseOnlyInitializedMappings;

    element_no_t elementNo = 0;
    int ghostMeshNo = 0;
    const int D = FunctionSpaceTargetType::dim();
    std::array<double,D> xi;

    // create the object that guesses the next target element no where the source dof should be found
    TargetElementNoEstimator<FunctionSpaceSourceType,FunctionSpaceTargetType> targetElementNoEstimator(functionSpaceSource_, functionSpaceTarget_);
    targetMappingInfo_.resize(nDofsLocalSource);

    std::vector<bool> targetDofIsMappedTo(nDofsLocalTarget, false);   //< for every target dof if it will get a value from any source dof

    if (VLOG_IS_ON(1))
    {
      VLOG(1) << "create mapping " << functionSpaceSource->meshName() << " -> " << functionSpaceTarget->meshName();
      VLOG(1) << "source geometry: " << functionSpaceSource->geometryField();
      VLOG(1) << "target geometry: " << functionSpaceTarget->geometryField();
      VLOG(1) << "target meshPartition: " << *functionSpaceTarget->meshPartition();
    }

    //VLOG(1) << "geometryField: " << functionSpaceTarget->geometryField();

    // xiTolerance is the tolerance for the element-local coordinate xi \in [0,1]^D, up to where a point is still considered and checked if it is inside the element.
    // A point is also inside the element for xi < 0 but xi > -xiTolerance or, analogously, xi > 1 but xi < 1+xiTolerance, if the respective element is the best fit,
    // i.e. if the point is outside of the actual mesh, then it is assumed (for the mapping) that it is inside the nearest element.
    // if xi tolerance was not set, set to default value
    if (xiTolerance <= 0)
      xiTolerance = 1e-2;
      
    bool startSearchInCurrentElement = true;    // start in element 0, maybe this is already the first element (it is if both meshes are completely aligned)
    int nSourceDofsOutsideTargetMesh = 0;
    double residual;
    bool searchedAllElements = false;
    int nTimesSearchedAllElements = 0;

    // visualization for 1D-1D: s=source, t=target, source dim <= target dim
    // t--s--------t-----s-----t

    // loop over all local dofs of the source functionSpace
    for (dof_no_t sourceDofNoLocal = 0; sourceDofNoLocal != nDofsLocalSource; sourceDofNoLocal++)
    {
      // determine information how to map a source value to the target mesh
      targetDof_t targetMappingInfo;
      targetMappingInfo.targetElements.resize(1);

      // Predict the next elementNo that will contain the sourceDof by a heuristic. If it is right, the search for the element is omitted, if it is wrong,
      // there will be a search first among the neighbouring elements,
      // if this was not successful (should not happen), then there is a search among all local elements.
      // In the special case where two 3D meshes are mapped onto each other and one mesh is a superset of the other, update the elementNo to the exact no.
      targetElementNoEstimator.estimateElementNo(sourceDofNoLocal, elementNo);

      // get node position of the source dof
      //dof_no_t sourceDofNoGlobal = functionSpaceTarget->meshPartition()->getDofNoGlobalPetsc(sourceDofNoLocal);
      Vec3 position = functionSpaceSource->getGeometry(sourceDofNoLocal);

      // find element no in the target mesh where the position is
      if (functionSpaceTarget->findPosition(position, elementNo, ghostMeshNo, xi, startSearchInCurrentElement, residual, searchedAllElements, xiTolerance))
      {
        targetMappingInfo.mapThisDof = true;

        // if there was a search among all elements, output a debugging message about this event
        if (searchedAllElements)
        {
          nTimesSearchedAllElements++;

          if (enableWarnings)
          {
            LOG(DEBUG) << "searched all elements to find elementNo: " << elementNo << " with xi: " << xi << ", xiTolerance: " << xiTolerance << ", residual: " << residual;
          }
        }

        VLOG(1) << "found at xi=" << xi << ", elementNo: " << elementNo << ", xiTolerance=" << xiTolerance << ", searchedAllElements: " << searchedAllElements << ", residual: " << residual;
      }
      else
      {
        if (enableWarnings)
        {
          LOG(INFO) << "In mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
            << functionSpaceTarget->meshName() << "\", source dof local " << sourceDofNoLocal
            << " of mesh \"" << functionSpaceSource->meshName() << "\" at position " << position << " is outside of target mesh \""
            << functionSpaceTarget->meshName() << "\" with tolerance " << xiTolerance << ". Try increasing parameter \"xiTolerance\".";
          LOG(INFO) << "position: " << position << ", startSearchInCurrentElement: " << startSearchInCurrentElement << ", xi: " << xi << ", residual: " << residual
            << ", elementNo: " << elementNo << ", searchedAllElements: " << searchedAllElements;
        }

        nSourceDofsOutsideTargetMesh++;
        targetMappingInfo.mapThisDof = false;
      }

      // store element no
      targetMappingInfo.targetElements[0].elementNoLocal = elementNo;

      std::array<dof_no_t,FunctionSpaceTargetType::nDofsPerElement()> targetDofNos = functionSpaceTarget->getElementDofNosLocal(elementNo);

      // determine factors how to distribute the source value to the dofs of the target element

      // note: geometry value = sum over dofs of geometryValue_dof * phi_dof(xi)
      for (int targetDofIndex = 0; targetDofIndex < nDofsPerTargetElement; targetDofIndex++)
      {
        double phiContribution;

        // for quadratic elements, treat as consisting of linear elements, this is disabled because it gives worse quality than the direct quadratic contributions
        if (false && std::is_same<typename FunctionSpaceTargetType::BasisFunction,typename BasisFunction::LagrangeOfOrder<2>>::value)
        {
          bool sourceDofHasContributionToTargetDof = true;
          phiContribution = quadraticElementComputePhiContribution(xi, targetDofIndex, sourceDofHasContributionToTargetDof);

          if (!sourceDofHasContributionToTargetDof)
            continue;
        }
        else
        {
          // for linear elements
          phiContribution = functionSpaceTarget->phi(targetDofIndex, xi);
        }

        // if phi is close to zero, set to 1e-14, this is practically zero, but it is still possible to divide by it in case the dof does not get any other contribution
        if (fabs(phiContribution) < 1e-14)
        {
          if (phiContribution >= 0)
          {
            phiContribution = 1e-14;
          }
          else
          {
            phiContribution = -1e-14;
          }
        }
        else
        {
          dof_no_t targetDofNoLocal = targetDofNos[targetDofIndex];

          // if this dof is local, store information that this target dof will get a value in the mapping,
          // i.e. there is a source dof that influences the mapped value of the target dof
          if (targetDofNoLocal < nDofsLocalTarget)
          {
            targetDofIsMappedTo[targetDofNoLocal] = true;
          }
        }

        targetMappingInfo.targetElements[0].scalingFactors[targetDofIndex] = phiContribution;
      }

      // debugging output about how interpolation is done, only in debug mode
#ifndef NDEBUG
      std::array<Vec3,FunctionSpaceTargetType::nDofsPerElement()> targetPositions;
      functionSpaceTarget->geometryField().getElementValues(targetMappingInfo.targetElements[0].elementNoLocal, targetPositions);
      Vec3 computedSourcePosition{0,0,0};

      for (int i = 0; i < FunctionSpaceTargetType::nDofsPerElement(); i++)
        computedSourcePosition += targetPositions[i] * targetMappingInfo.targetElements[0].scalingFactors[i];

      // print error for mismatch
      if (fabs(computedSourcePosition[0] - position[0]) > 1e-1 || fabs(computedSourcePosition[1] - position[1]) > 1e-1 || fabs(computedSourcePosition[2] - position[2]) > 1e-1)
        LOG(ERROR) << "mismatch " << computedSourcePosition << " != " << position << ", " << ", xi: " << xi;

      LOG(DEBUG) << functionSpaceSource->meshName() << " dof " << sourceDofNoLocal << " value = " << position << " = " << targetPositions << " * " << targetMappingInfo.targetElements[0].scalingFactors
        << ", interpolation in element " << targetMappingInfo.targetElements[0].elementNoLocal << " of " << functionSpaceTarget->meshName();
#endif


      targetMappingInfo_[sourceDofNoLocal] = targetMappingInfo;

      if (VLOG_IS_ON(2))
      {
        double scalingFactorsSum = 0;
        for (int targetDofIndex = 0; targetDofIndex < nDofsPerTargetElement; targetDofIndex++)
        {
          scalingFactorsSum += targetMappingInfo_[sourceDofNoLocal].targetElements[0].scalingFactors[targetDofIndex];
        }
        VLOG(2) << "  source dof local " << sourceDofNoLocal << ", pos: " << position << ", xi: " << xi
          << ", element no: " << targetMappingInfo.targetElements[0].elementNoLocal << ", scaling factors: " << targetMappingInfo_[sourceDofNoLocal].targetElements[0].scalingFactors
          << ", sum: " << scalingFactorsSum;
      }
    }

    int nTargetDofsNotMapped = 0;
    int nTimesSearchedAllElementsForFix = 0;
    int nTargetDofNosLocaNotFixed = 0;

    // find target dofs that do not appear in any targetMappingInfo and therefore will so far not receive any value when mapping from source to target
    fixUnmappedDofs(functionSpaceSource, functionSpaceTarget, xiTolerance, compositeUseOnlyInitializedMappings, isEnabledFixUnmappedDofs, targetDofIsMappedTo,
                    nTargetDofsNotMapped, nTimesSearchedAllElementsForFix, nTargetDofNosLocaNotFixed);

    Control::PerformanceMeasurement::stop("durationComputeMappingBetweenMeshes");

    if (nSourceDofsOutsideTargetMesh > 0)
    {
      LOG(INFO) << "Successfully initialized mapping between meshes \"" << functionSpaceSource->meshName() << "\" and \""
        << functionSpaceTarget->meshName() << "\", " << nSourceDofsOutsideTargetMesh << "/" << nDofsLocalSource << " source dofs are outside the target mesh. "
        << (!enableWarnings ? "\"enableWarnings: False\"" : "\"enableWarnings: True\"")
        << " \"xiTolerance\": " << xiTolerance << ", total duration of all mappings: " << Control::PerformanceMeasurement::getDuration("durationComputeMappingBetweenMeshes") << " s";
    }

    // add statistics to log

    // add log message, to be included in the log file
    std::stringstream logMessage;
    logMessage << "  Statistics: " << nSourceDofsOutsideTargetMesh << "/" << nDofsLocalSource << " source dofs are outside the target mesh,\n";
    if (nTimesSearchedAllElementsForFix != 0)
    {
      logMessage << "              " << nTargetDofNosLocaNotFixed << "/" << nDofsLocalTarget << " target dofs are outside the source mesh, "
        << "thereof " << nTargetDofsNotMapped-nTargetDofNosLocaNotFixed << " previously unmapped target dofs were added,\n"
        << "              iterated " << nTimesSearchedAllElements << " times over target mesh for mapping, "
        << nTimesSearchedAllElementsForFix << " times for fixing unmapped target dofs,\n"
        << "              \"xiTolerance\": " << xiTolerance << " (increase this value to reduce the number of (costly) iterations over the whole mesh, however increasing potentially leads to more elements being checked which takes longer).\n";
    }
    else
    {
      logMessage << "              " << nTargetDofsNotMapped << "/" << nDofsLocalTarget << " target dofs are outside the source mesh (fixing this is disabled),\n"
        << "              iterated " << nTimesSearchedAllElements << " times over target mesh,\n"
        << "              \"xiTolerance\": " << xiTolerance << " (increase this value to reduce the number of (costly) iterations over the whole mesh, however increasing potentially leads to more elements being checked which takes longer).\n";
    }
    logMessage << "              Total duration of all mappings so far: " << Control::PerformanceMeasurement::getDuration("durationComputeMappingBetweenMeshes") << " s.";
    DihuContext::mappingBetweenMeshesManager()->addLogMessage(logMessage.str());

  }  // if not composite
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
double MappingBetweenMeshesConstruct<FunctionSpaceSourceType, FunctionSpaceTargetType>::
quadraticElementComputePhiContribution(std::array<double,FunctionSpaceTargetType::dim()> xi,
                                       int targetDofIndex, bool &sourceDofHasContributionToTargetDof)
{
  sourceDofHasContributionToTargetDof = true;

  // split the quadratic element with 3^D nodes into 2^D sub elements with 2^D nodes each
  const int D = FunctionSpaceTargetType::dim();
  std::array<double,D> xiSubElement;    //< xi value of the target element
  int targetDofIndexSubElement = 0;     //< dof index of the target sub element

  // loop over coordinate directions
  for (int coordinateDirectionNo = 0; coordinateDirectionNo < D; coordinateDirectionNo++)
  {
    int targetDofIndexCoordinateDirection = 0;
    switch (coordinateDirectionNo)
    {
    case 0:
      targetDofIndexCoordinateDirection = targetDofIndex % 3;
      break;

    case 1:
      targetDofIndexCoordinateDirection = int(targetDofIndex % 9 / 3);
      break;

    case 2:
      targetDofIndexCoordinateDirection = int(targetDofIndex / 9);
      break;
    }

    // if the source dof is in the first half of the quadratic element along the current coordinate direction
    if (xi[coordinateDirectionNo] < 0.5)
    {
      // if the target dof is in the second half of the quadratic element along the current coordinate direction
      if (targetDofIndexCoordinateDirection == 2)
      {
        sourceDofHasContributionToTargetDof = false;
        break;
      }

      targetDofIndexSubElement += pow(2,coordinateDirectionNo) * targetDofIndexCoordinateDirection;

      // compute the xi value of the sub element
      xiSubElement[coordinateDirectionNo] = xi[coordinateDirectionNo] * 2;
    }
    else
    {
      // if the source dof is in the second half and the targetDof is in the first half, we do not have a contributino
      if (targetDofIndexCoordinateDirection == 0)
      {
        sourceDofHasContributionToTargetDof = false;
        break;
      }

      targetDofIndexSubElement += pow(2,coordinateDirectionNo) * (targetDofIndexCoordinateDirection - 1);

      // compute the xi value of the sub element
      xiSubElement[coordinateDirectionNo] = (xi[coordinateDirectionNo]-0.5) * 2;
    }
  }

  double phiContribution = FunctionSpace::FunctionSpaceFunction<typename FunctionSpaceTargetType::Mesh, BasisFunction::LagrangeOfOrder<1>>::phi(targetDofIndexSubElement, xiSubElement);

  if (sourceDofHasContributionToTargetDof)
    LOG(INFO) << "dof: " << targetDofIndex << " -> " << targetDofIndexSubElement << ", xi: " << xi << " -> " << xiSubElement << " phiContribution: " << phiContribution;

  return phiContribution;
}

template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
void MappingBetweenMeshesConstruct<FunctionSpaceSourceType, FunctionSpaceTargetType>::
fixUnmappedDofs(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget,
                double xiTolerance, bool compositeUseOnlyInitializedMappings, bool isEnabledFixUnmappedDofs, const std::vector<bool> &targetDofIsMappedTo,
                int &nTargetDofsNotMapped, int &nTimesSearchedAllElements, int &nTargetDofNosLocaNotFixed)
{
  const dof_no_t nDofsLocalTarget = functionSpaceTarget->nDofsLocalWithoutGhosts();
  const int nDofsPerTargetElement = FunctionSpaceTargetType::nDofsPerElement();

  // count number of target dofs that have no source dof mapped
  nTargetDofsNotMapped = 0;
  nTimesSearchedAllElements = 0;

  for (dof_no_t targetDofNoLocal = 0; targetDofNoLocal < nDofsLocalTarget; targetDofNoLocal++)
  {
    if (!targetDofIsMappedTo[targetDofNoLocal])
    {
      nTargetDofsNotMapped++;
    }
  }

  // check if the algorithm should be performed here at all
  if (FunctionSpaceSourceType::dim() >= FunctionSpaceTargetType::dim() && nTargetDofsNotMapped > 0 && compositeUseOnlyInitializedMappings)
  {
    LOG(DEBUG) << "mapping \"" << functionSpaceSource->meshName() << "\" -> \""
      << functionSpaceTarget->meshName() << "\", source FunctionSpace dim: "
      << FunctionSpaceSourceType::dim() << " >= target FunctionSpace dim: " << FunctionSpaceTargetType::dim() << ", "
      << nTargetDofsNotMapped << " target dofs of " << nDofsLocalTarget << " have no source dofs that would contribute values. "
      << "But option \"compositeUseOnlyInitializedMappings\" is set to True, therefore not fixing the missing target dofs (as they might be on a different submesh)";

    // add log message, to be included in the log file
    std::stringstream logMessage;
    logMessage << "  " << nTargetDofsNotMapped << " target dofs of " << nDofsLocalTarget << " have no source dofs that would contribute values. \n"
      << "  But option \"compositeUseOnlyInitializedMappings\" is set to True, therefore not fixing the missing target dofs (as they might be on a different submesh)";

    DihuContext::mappingBetweenMeshesManager()->addLogMessage(logMessage.str());
    return;
  }
  else if (!isEnabledFixUnmappedDofs)
  {
    return;
  }

  if (FunctionSpaceSourceType::dim() >= FunctionSpaceTargetType::dim() && nTargetDofsNotMapped > 0 && !compositeUseOnlyInitializedMappings)
  {
    // do the algorithm

    // debugging output
    LOG(DEBUG) << "mapping \"" << functionSpaceSource->meshName() << "\" -> \""
      << functionSpaceTarget->meshName() << "\":" << nTargetDofsNotMapped << " target dofs of " << nDofsLocalTarget
      << " have no source dofs "
      << "that would contribute values. Source FunctionSpace dim: "
      << FunctionSpaceSourceType::dim() << " >= target FunctionSpace dim: " << FunctionSpaceTargetType::dim() << ". Now fixing.";

    std::set<dof_no_t> targetDofNoLocalNotFixed;    // collect all dofs that are still not fixed
    nTimesSearchedAllElements = 1;

    // loop over all elements in the target function space
    for (element_no_t targetElementNoLocal = 0; targetElementNoLocal < functionSpaceTarget->nElementsLocal(); targetElementNoLocal++)
    {
      std::array<dof_no_t,FunctionSpaceTargetType::nDofsPerElement()> targetDofNos = functionSpaceTarget->getElementDofNosLocal(targetElementNoLocal);

      // loop over dofs of target element
      for (int targetDofIndex = 0; targetDofIndex < nDofsPerTargetElement; targetDofIndex++)
      {
        dof_no_t targetDofNoLocal = targetDofNos[targetDofIndex];

        // if this is a ghost dof, do not handle it
        if (targetDofNoLocal >= nDofsLocalTarget)
          continue;

        Vec3 position = functionSpaceTarget->getGeometry(targetDofNoLocal);
        LOG(DEBUG) << " e" << targetElementNoLocal << " i" << targetDofIndex << ", targetDofIsMappedTo[" << targetDofNoLocal << "]: " << targetDofIsMappedTo[targetDofNoLocal] << ", position: " << position;

        // if target dof is not being mapped to by any source dof, simply initiate interpolation of the source mesh to this dof
        if (!targetDofIsMappedTo[targetDofNoLocal])
        {
          // get one element of the target dof

          // find element and xi position in source mesh where target dof is located
          Vec3 position = functionSpaceTarget->getGeometry(targetDofNoLocal);
          element_no_t sourceElementNo = 0;
          bool startSearchInCurrentElement = false;
          std::array<double,FunctionSpaceSourceType::dim()> xiSource;
          int ghostMeshNo = 0;
          double residual;
          bool searchedAllElements = false;

          //LOG(DEBUG) << "target (el." << targetElementNoLocal << ",index" << targetDofIndex << ") dof " << targetDofNoLocal << ", position: " << position
          //  << " is not mapped, now find element in source function space";
          if (functionSpaceSource->findPosition(position, sourceElementNo, ghostMeshNo, xiSource, startSearchInCurrentElement, residual, searchedAllElements, xiTolerance))
          {
            if (searchedAllElements)
              nTimesSearchedAllElements++;

            // get dofs of this source element
            std::array<dof_no_t,FunctionSpaceSourceType::nDofsPerElement()> sourceDofNos = functionSpaceSource->getElementDofNosLocal(sourceElementNo);

            LOG(DEBUG) << "at position " << position << " found source element " << sourceElementNo << ", xi " << xiSource << ", with dofs " << sourceDofNos;

            // loop over all the source dofs that will contribute to the value of the target dof
            for (int sourceDofIndex = 0; sourceDofIndex != FunctionSpaceSourceType::nDofsPerElement(); sourceDofIndex++)
            {
              dof_no_t sourceDofNoLocal = sourceDofNos[sourceDofIndex];

              // create new entry for the targetMappingInfo_[sourceDofNoLocal]
              typename targetDof_t::element_t targetElement;

              // set element no
              targetElement.elementNoLocal = targetElementNoLocal;

              // set scaling factors
              for (int i = 0; i < FunctionSpaceTargetType::nDofsPerElement(); i++)
              {
                targetElement.scalingFactors[i] = 0;
              }

              double phiContribution = functionSpaceSource->phi(sourceDofIndex, xiSource);
              targetElement.scalingFactors[targetDofIndex] = phiContribution;

              try
              {
                targetMappingInfo_[sourceDofNoLocal].targetElements.push_back(targetElement);
              }
              catch (...)
              {
                LOG(ERROR) << "Could allocate memory while creation of mapping \"" << functionSpaceSource->meshName() << "\" -> \""
                  << functionSpaceTarget->meshName() << "\":" << nTargetDofsNotMapped << " target dofs of " << nDofsLocalTarget
                  << " have no source dofs that would contribute values. Source FunctionSpace dim: "
                  << FunctionSpaceSourceType::dim() << " >= target FunctionSpace dim: " << FunctionSpaceTargetType::dim();
                break;
              }

              //LOG(DEBUG) << "add scaling Factor " << phiContribution << " at targetDofIndex " << targetDofIndex << " of targetELement " << targetElementNoLocal
              //  << ", now, source dof " << sourceDofNoLocal << " has " << targetMappingInfo_[sourceDofNoLocal].targetElements.size() << " target elements.";

            }

            // now the target dof is fixed
            //targetDofIsMappedTo[targetDofNoLocal] = true;
          }
          else
          {
            LOG(DEBUG) << "Could not find element of source dof for position " << position;
            targetDofNoLocalNotFixed.insert(targetDofNoLocal);
          }
        }
      }
    }
    LOG(DEBUG) << "after fixing target dofs by source mesh interpolation, " << targetDofNoLocalNotFixed.size()
      << " remaining targetDofNoLocalNotFixed: " << targetDofNoLocalNotFixed;

    // add log message, to be included in the log file
    std::stringstream logMessage;
    logMessage << "  " << nTargetDofsNotMapped << " target dofs of " << nDofsLocalTarget << " had no source dofs that would contribute values."
      << "Option \"fixUnmappedDofs\" is set to True. After source mesh interpolation, " << targetDofNoLocalNotFixed.size() << " target dofs are still unmapped. "
      << "nTimesSearchedAllElements: " << nTimesSearchedAllElements;

    DihuContext::mappingBetweenMeshesManager()->addLogMessage(logMessage.str());

    nTargetDofNosLocaNotFixed = targetDofNoLocalNotFixed.size();
  }
  else
  {
    LOG(DEBUG) << nTargetDofsNotMapped << " target dofs have no source dofs that would contribute values. Source FunctionSpace dim: "
      << FunctionSpaceSourceType::dim() << ", target FunctionSpace dim: " << FunctionSpaceTargetType::dim() << ".";
  }
}

//! get access to the internal targetMappingInfo_ variable
template<typename FunctionSpaceTargetType, typename FunctionSpaceSourceType>
const std::vector<typename MappingBetweenMeshesConstruct<FunctionSpaceTargetType, FunctionSpaceSourceType>::targetDof_t> &
MappingBetweenMeshesConstruct<FunctionSpaceTargetType, FunctionSpaceSourceType>::
targetMappingInfo() const
{
  return targetMappingInfo_;
}

}  // namespace

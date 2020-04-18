#include "mesh/mapping_between_meshes/mapping/02_composite.h"

#include "control/diagnostic_tool/performance_measurement.h"

#include "utility/vector_operators.h"
#include "control/dihu_context.h"
#include "mesh/mesh_manager/mesh_manager.h"
#include "mesh/mapping_between_meshes/manager/01_manager.h"

namespace MappingBetweenMeshes
{

template<int D, typename BasisFunctionType, typename FunctionSpaceTargetType>
MappingBetweenMeshes<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>, FunctionSpaceTargetType>::
MappingBetweenMeshes(std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>> functionSpaceSource,
                              std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget,
                              double xiTolerance, bool enableWarnings, bool compositeUseOnlyInitializedMappings) :
  MappingBetweenMeshesImplementation<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>, FunctionSpaceTargetType>(
    functionSpaceSource, functionSpaceTarget, xiTolerance, enableWarnings, compositeUseOnlyInitializedMappings)
{
  if (compositeUseOnlyInitializedMappings)
  {
    LOG(DEBUG) << "Composite mesh \"" << this->functionSpaceSource_->meshName() << "\", mapping to \"" << this->functionSpaceTarget_->meshName() << "\".";
    typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> SubFunctionSpaceType;

    // get the sub function spaces of the composite function spaces
    const std::vector<std::shared_ptr<SubFunctionSpaceType>> &sourceSubFunctionSpaces
      = this->functionSpaceSource_->subFunctionSpaces();

    // make types available    
    typedef typename MappingBetweenMeshesImplementation<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>, FunctionSpaceTargetType>::targetDof_t targetDof_t;
    typedef typename MappingBetweenMeshesImplementation<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>, FunctionSpaceTargetType>::targetDof_t::element_t element_t;

    targetDof_t newTargetDof;
    newTargetDof.mapThisDof = false;
    this->targetMappingInfo_.resize(this->functionSpaceSource_->nDofsLocalWithoutGhosts(), newTargetDof);

    // loop over submeshes
    for (int subMeshNo = 0; subMeshNo < sourceSubFunctionSpaces.size(); subMeshNo++)
    {
      // get current sub meshes
      std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> subFunctionSpace = sourceSubFunctionSpaces[subMeshNo];

      // if the mapping exists
      if (DihuContext::mappingBetweenMeshesManager()->hasMappingBetweenMeshes<SubFunctionSpaceType,FunctionSpaceTargetType>(subFunctionSpace, this->functionSpaceTarget_))
      {
        // get the mapping
        std::shared_ptr<MappingBetweenMeshes<SubFunctionSpaceType,FunctionSpaceTargetType>> mapping 
          = DihuContext::mappingBetweenMeshesManager()->mappingBetweenMeshes<SubFunctionSpaceType,FunctionSpaceTargetType>(subFunctionSpace, this->functionSpaceTarget_);

        LOG(DEBUG) << "add information from mapping \"" << subFunctionSpace->meshName() << "\"->\"" << this->functionSpaceTarget_->meshName() << "\" to mapping "
          << "\"" << this->functionSpaceSource_->meshName() << "\"->\"" << this->functionSpaceTarget_->meshName() << "\".";

        if (!mapping)
        {
          LOG(FATAL) << "mapping is not set";
        }

        // add mapping to own mapping
        dof_no_t sourceDofNoLocalOnSubMesh = 0;
        for (typename std::vector<typename MappingBetweenMeshes<SubFunctionSpaceType,FunctionSpaceTargetType>::targetDof_t>::const_iterator targetDofIter = mapping->targetMappingInfo().begin();
             targetDofIter != mapping->targetMappingInfo().end(); targetDofIter++, sourceDofNoLocalOnSubMesh++)
        {
          // create a new target dof object that will be added to the composite mapping
          targetDof_t newTargetDof;
          newTargetDof.mapThisDof = targetDofIter->mapThisDof;

          // loop over target elements of the current dof in the mapping
          for (auto &targetElement : targetDofIter->targetElements)
          {
            // fill a new targetElement for the composite mapping with the information from the mapping of the submesh
            element_t newTargetElement;
            newTargetElement.elementNoLocal = targetElement.elementNoLocal;
            newTargetElement.scalingFactors = targetElement.scalingFactors;
            newTargetDof.targetElements.push_back(newTargetElement);

            if (newTargetElement.elementNoLocal >= this->functionSpaceTarget_->nElementsLocal())
            {
              LOG(FATAL) << "in creating composite mapping \"" << this->functionSpaceSource_->meshName() << "\" -> \"" << this->functionSpaceTarget_->meshName() 
                << "\", mapping  \"" << subFunctionSpace->meshName() << "\"->\"" << this->functionSpaceTarget_->meshName() << "\", targetElement " << newTargetElement.elementNoLocal
                << " is out of range, number of local elements in \"" << this->functionSpaceTarget_->meshName() << "\": " << this->functionSpaceTarget_->nElementsLocal();
            }
          }

          // convert the source dof no in the sub mesh numbering to the source dof no in the composite numbering
          const int nDofsPerNodeSource = functionSpaceSource->nDofsPerNode();
          int dofIndex = sourceDofNoLocalOnSubMesh % nDofsPerNodeSource;
          node_no_t sourceNodeNoLocalOnSubMesh = sourceDofNoLocalOnSubMesh / nDofsPerNodeSource;

          // from the submesh no and the local node no in the submesh numbering get the local node no in the composite numbering
          bool nodeIsSharedAndRemovedInCurrentMesh;
          node_no_t sourceNodeNoLocalComposite = this->functionSpaceSource_->meshPartition()->getNodeNoLocalFromSubmesh(subMeshNo, sourceNodeNoLocalOnSubMesh, nodeIsSharedAndRemovedInCurrentMesh);
          dof_no_t sourceDofNoLocalComposite = sourceNodeNoLocalComposite * nDofsPerNodeSource + dofIndex;

          // add information to own mapping
          this->targetMappingInfo_[sourceDofNoLocalComposite] = newTargetDof;

          if (this->targetMappingInfo_.back().targetElements.empty())
            LOG(DEBUG) << "source dof " << sourceDofNoLocalComposite << " (submesh " << subMeshNo << " local dof no " << sourceDofNoLocalOnSubMesh << "), no target elements";
          else
            LOG(DEBUG) << "source dof " << sourceDofNoLocalComposite << " (submesh " << subMeshNo << " local dof no " << sourceDofNoLocalOnSubMesh << ")"
               << ", element no local " << this->targetMappingInfo_[sourceDofNoLocalComposite].targetElements[0].elementNoLocal
               << ", scaling factors " << this->targetMappingInfo_[sourceDofNoLocalComposite].targetElements[0].scalingFactors;
        }
      }
      else 
      {
        LOG(DEBUG) << "No mapping \"" << subFunctionSpace->meshName() << "\"->\"" << this->functionSpaceTarget_->meshName() << "\" exists.";
      }
    }  // loop over submeshes

#ifndef NDEBUG

    // debugging output

    // collect statistics
    int nEntriesWithNoTargetElements = 0;
    int nEntriesWithOneTargetElement = 0;
    int nEntriesWithMultipleTargetElements = 0;
    int nUnmappedEntries = 0;

    for (dof_no_t sourceDofNoLocal = 0; sourceDofNoLocal < this->targetMappingInfo_.size(); sourceDofNoLocal++)
    {
      if (this->targetMappingInfo_[sourceDofNoLocal].targetElements.size() == 0)
        nEntriesWithNoTargetElements++;
    
      else if (this->targetMappingInfo_[sourceDofNoLocal].targetElements.size() == 1)
        nEntriesWithOneTargetElement++;

      else 
        nEntriesWithMultipleTargetElements++;

      if (!this->targetMappingInfo_[sourceDofNoLocal].mapThisDof)
        nUnmappedEntries++;

      LOG(DEBUG) << this->targetMappingInfo_.size() << " target mapping dofs, " << nUnmappedEntries << " unmapped, " << nEntriesWithNoTargetElements 
        << " with no target elements, " << nEntriesWithOneTargetElement << " with 1 target element, " << nEntriesWithMultipleTargetElements << " with multiple target elements.";
    }

#endif

  } // if compositeUseOnlyInitializedMappings
}

}  // namespace

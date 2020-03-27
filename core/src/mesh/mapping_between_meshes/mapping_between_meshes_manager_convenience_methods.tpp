#include "mesh/mapping_between_meshes/mapping_between_meshes_manager.h"

namespace Mesh
{

/** helper classes */
template<typename FieldVariableSourceType, typename FieldVariableTargetType, typename Dummy=void>
struct ExtractComponentShared
{
  // helper function, does nothing
  static void call(int componentNoSource, std::shared_ptr<FieldVariableSourceType> fieldVariableSource, std::shared_ptr<FieldVariableTargetType> fieldVariableTarget)
  {
    LOG(FATAL) << "This method should not be called, something went wrong! Target field variable has "
      << FieldVariableTargetType::nComponents() << " != 1 components or dimension mismatch: "
      << FieldVariableSourceType::FunctionSpace::dim() << " != "<< FieldVariableTargetType::FunctionSpace::dim();
  }
};

template<typename FieldVariableSourceType, typename FieldVariableTargetType>
struct ExtractComponentShared<FieldVariableSourceType,FieldVariableTargetType,
  typename std::enable_if<FieldVariableTargetType::nComponents() == 1 && FieldVariableSourceType::FunctionSpace::dim() == FieldVariableTargetType::FunctionSpace::dim() 
     && std::is_same<typename FieldVariableSourceType::FunctionSpace::BasisFunction,typename FieldVariableTargetType::FunctionSpace::BasisFunction>::value,void>::type
>
{
  // extract the component from the field variable
  static void call(int componentNoSource, std::shared_ptr<FieldVariableSourceType> fieldVariableSource, std::shared_ptr<FieldVariableTargetType> fieldVariableTarget)
  {
    fieldVariableSource->extractComponentShared(componentNoSource, fieldVariableTarget);
  }
};

template<typename FieldVariableSourceType, typename FieldVariableTargetType, typename Dummy=void>
struct ExtractComponentCopy
{
  // helper function, does nothing
  static void call(int componentNoSource, std::shared_ptr<FieldVariableSourceType> fieldVariableSource, std::shared_ptr<FieldVariableTargetType> fieldVariableTarget)
  {
    LOG(FATAL) << "This method should not be called, something went wrong! Target field variable has "
      << FieldVariableTargetType::nComponents() << " != 1 components or dimension mismatch: "
      << FieldVariableSourceType::FunctionSpace::dim() << " != "<< FieldVariableTargetType::FunctionSpace::dim();
  }
};

template<typename FieldVariableSourceType, typename FieldVariableTargetType>
struct ExtractComponentCopy<FieldVariableSourceType,FieldVariableTargetType,
  typename std::enable_if<FieldVariableTargetType::nComponents() == 1 && FieldVariableSourceType::FunctionSpace::dim() == FieldVariableTargetType::FunctionSpace::dim() 
     && std::is_same<typename FieldVariableSourceType::FunctionSpace::BasisFunction,typename FieldVariableTargetType::FunctionSpace::BasisFunction>::value,void>::type
>
{
  // extract the component from the field variable
  static void call(int componentNoSource, std::shared_ptr<FieldVariableSourceType> fieldVariableSource, std::shared_ptr<FieldVariableTargetType> fieldVariableTarget)
  {
    fieldVariableSource->extractComponentCopy(componentNoSource, fieldVariableTarget);
  }
};


template<typename FieldVariableSourceType, typename FieldVariableTargetType, typename Dummy=void>
struct RestoreExtractedComponent
{
  static void call(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget)
  {
    LOG(FATAL) << "This method should not be called, something went wrong! Dimension of field variables mismatches, "
      << FieldVariableTargetType::FunctionSpace::dim() << "!=" << FieldVariableSourceType::FunctionSpace::dim();
  }
};

template<typename FieldVariableSourceType, typename FieldVariableTargetType>
struct RestoreExtractedComponent<FieldVariableSourceType,FieldVariableTargetType,
  typename std::enable_if<FieldVariableTargetType::FunctionSpace::dim() == FieldVariableSourceType::FunctionSpace::dim() 
     && std::is_same<typename FieldVariableSourceType::FunctionSpace::BasisFunction,typename FieldVariableTargetType::FunctionSpace::BasisFunction>::value,void>::type
>
{
  static void call(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget)
  {
    fieldVariableTarget->restoreExtractedComponent(fieldVariableSource->partitionedPetscVec(), componentNoTarget);
  };
};


//! Simplified methods

//! map between meshes of any dimensionality,
//! if the field variables share the same function space, do no mapping at all, but copy the values
template<typename FieldVariableSourceType, typename FieldVariableTargetType>
void MappingBetweenMeshesManager::
map(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
    std::shared_ptr<FieldVariableTargetType> &fieldVariableTarget, int componentNoTarget, bool avoidCopyIfPossible)
{
  VLOG(1) << "map " << fieldVariableSource->name() << "." << componentNoSource << " (dim " << FieldVariableSourceType::FunctionSpace::dim() << ", "
    << FieldVariableSourceType::nComponents() << " components total)"
    << " -> " << fieldVariableTarget->name() << "." << componentNoTarget << " (dim " << FieldVariableTargetType::FunctionSpace::dim() << ", "
    << FieldVariableTargetType::nComponents() << " components total)"
    << "), avoidCopyIfPossible: " << avoidCopyIfPossible;

  // assert that both or none of the componentNos are -1
  assert((componentNoSource == -1) == (componentNoTarget == -1));

  // At first, test if the function space is the same and the values can simply be copied or extracted.
  // if the dimensionality of source and target function space is the same
  if (FieldVariableSourceType::FunctionSpace::dim() == FieldVariableTargetType::FunctionSpace::dim())
  {
    // check if the function space is also the same (i.e., same mesh)
    if (fieldVariableSource->functionSpace()->meshName() == fieldVariableTarget->functionSpace()->meshName())
    {
      VLOG(1) << "mesh is the same: " << fieldVariableSource->functionSpace()->meshName() << " -> " << fieldVariableTarget->functionSpace()->meshName();
      VLOG(1) << "target representation: " << fieldVariableTarget->partitionedPetscVec()->getCurrentRepresentationString();

      // if representation of fieldVariableTarget is invalid, this means that it has been extracted to another field variable
      if (fieldVariableTarget->partitionedPetscVec()->currentRepresentation() == Partition::values_representation_t::representationInvalid)
      {
        VLOG(1) << "call restore extracted component " << componentNoTarget;

        // transfer back, e.g. from finite elements back to cellml
        //fieldVariableTarget->restoreExtractedComponent(fieldVariableSource->partitionedPetscVec());
        RestoreExtractedComponent<FieldVariableSourceType,FieldVariableTargetType>::call(fieldVariableSource, fieldVariableTarget, componentNoTarget);
      }
      else
      {
        if (FieldVariableTargetType::nComponents() == 1)
        {
          if (avoidCopyIfPossible)
          {
            if (FieldVariableSourceType::nComponents() == 1)
            {
              // Here, both field variables are scalar and the mesh is the same and copy should be avoided.
              // This means that the target field variable pointer can be set to the source field variable pointer.
              if ((long long)(fieldVariableSource.get()) != (long long)(fieldVariableTarget.get()))
              {
                VLOG(1) << "set target pointer to source pointer";
                fieldVariableTarget = std::dynamic_pointer_cast<FieldVariableTargetType>(fieldVariableSource);
              }
              else
              {
                VLOG(1) << "Field variable is already the same";
              }
            }
            else if (fieldVariableSource->isExtractComponentSharedPossible(componentNoSource))
            {
              // fieldVariableTarget has only 1 component
              // The following retrieves the raw memory pointer from the Petsc vector in fieldVariableSource and reuses it for fieldVariableTarget
              // That means that fieldVariableSource cannot be used anymore, only after restoreExtractedComponent was called on fieldVariableSource.
              // This is done when mapping back from fieldVariableTarget to fieldVariableSource.
              if (componentNoSource == -1)
                componentNoSource = 0;

              VLOG(1) << "call extract component shared";

              //fieldVariableSource->extractComponentShared(componentNoSource, fieldVariableTarget);
              ExtractComponentShared<FieldVariableSourceType,FieldVariableTargetType>::call(componentNoSource, fieldVariableSource, fieldVariableTarget);

              VLOG(1) << "source representation: " << fieldVariableSource->partitionedPetscVec()->getCurrentRepresentationString();

              return;
            }
            else
            {
              // extractComponentShared is not possible, because the last component cannot be extracted by this
              if (componentNoSource == -1)
                componentNoSource = 0;

              VLOG(1) << "source component no " << componentNoSource << " of " << FieldVariableSourceType::nComponents() << ", call extract component copy";

              //fieldVariableSource->extractComponentCopy(componentNoSource, fieldVariableTarget);
              ExtractComponentCopy<FieldVariableSourceType,FieldVariableTargetType>::call(componentNoSource, fieldVariableSource, fieldVariableTarget);

              // if target field variable is "additionalFieldVariable", also set the name from the source field variable
              if (fieldVariableTarget->name().find("additionalFieldVariable") != std::string::npos)
              {
                fieldVariableTarget->setName(fieldVariableSource->name());
              }

              VLOG(1) << "source representation: " << fieldVariableSource->partitionedPetscVec()->getCurrentRepresentationString();

              return;
            }
          }
        }

        // fieldVariableTarget has > 1 components or avoidCopyIfPossible is false which means explicit copy
        if (componentNoSource == -1)
        {
          VLOG(1) << "copy all " << std::min(FieldVariableSourceType::nComponents(), FieldVariableTargetType::nComponents()) << " components";

          // Here, we copy the all components of fieldVariableSource to the corresponding components of fieldVariableTarget.
          PetscErrorCode ierr;
          for (int componentNo = 0; componentNo < std::min(FieldVariableSourceType::nComponents(), FieldVariableTargetType::nComponents()); componentNo++)
          {
            ierr = VecCopy(fieldVariableSource->valuesGlobal(componentNo), fieldVariableTarget->valuesGlobal(componentNo)); CHKERRV(ierr);
          }
        }
        else
        {
          VLOG(1) << "copy one component";

          // Here, we copy the given component of fieldVariableSource to the componentNoTarget of fieldVariableTarget.
          PetscErrorCode ierr;
          ierr = VecCopy(fieldVariableSource->valuesGlobal(componentNoSource), fieldVariableTarget->valuesGlobal(componentNoTarget)); CHKERRV(ierr);
        }
      }
      VLOG(1) << "source (" << fieldVariableSource << "): " << *fieldVariableSource;
      VLOG(1) << "target (" << fieldVariableTarget << "): " << *fieldVariableTarget;
      return;
    }
  }

  VLOG(1) << "map mesh " << fieldVariableSource->functionSpace()->meshName() << " -> " << fieldVariableTarget->functionSpace()->meshName();

  // The values could not be copied, call the appropriate mapping
  // if the dimensionality of source and target function space is the same
  if (FieldVariableSourceType::FunctionSpace::dim() > FieldVariableTargetType::FunctionSpace::dim())
  {
    VLOG(1) << "map high to low dimension";

    // this is the "inverse mapping" of a linear/bilinear/trilinear mapping
    mapHighToLowDimension(fieldVariableSource, componentNoSource, fieldVariableTarget, componentNoTarget);
  }
  else
  {
    VLOG(1) << "map low to high dimension";

    // this is the mapping that needs prepareMapping and finalizeMapping, it is the more natural one (linear/bilinear/trilinear)
    mapLowToHighDimension(fieldVariableSource, componentNoSource, fieldVariableTarget, componentNoTarget);
  }
}

//! prepare the mapping for meshes of any dimensionality, this can be called even if not needed
template<typename FieldVariableSourceType, typename FieldVariableTargetType>
void MappingBetweenMeshesManager::
prepareMapping(std::shared_ptr<FieldVariableSourceType> fieldVariableSource,
               std::shared_ptr<FieldVariableTargetType> fieldVariableTarget)
{
  // check if prepareMapping is needed
  if (FieldVariableSourceType::FunctionSpace::dim() == FieldVariableTargetType::FunctionSpace::dim())
  {
    // check if the function space is also the same (i.e., same mesh)
    if (fieldVariableSource->functionSpace()->meshName() == fieldVariableTarget->functionSpace()->meshName())
    {
      return;
    }
  }
  else if (FieldVariableSourceType::FunctionSpace::dim() > FieldVariableTargetType::FunctionSpace::dim())
  {
    return;
  }

  // prepareMapping is needed
  prepareMappingLowToHigh(fieldVariableTarget);
}

//! finalize the mapping for meshes of any dimensionality, this can be called even if not needed
template<typename FieldVariableSourceType, typename FieldVariableTargetType>
void MappingBetweenMeshesManager::
finalizeMapping(std::shared_ptr<FieldVariableSourceType> fieldVariableSource,
                std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget)
{
  // check if finalizeMapping is needed
  if (FieldVariableSourceType::FunctionSpace::dim() == FieldVariableTargetType::FunctionSpace::dim())
  {
    // check if the function space is also the same (i.e., same mesh)
    if (fieldVariableSource->functionSpace()->meshName() == fieldVariableTarget->functionSpace()->meshName())
    {
      return;
    }
  }
  else if (FieldVariableSourceType::FunctionSpace::dim() > FieldVariableTargetType::FunctionSpace::dim())
  {
    return;
  }

  // finalizeMapping is needed
  finalizeMappingLowToHigh(fieldVariableTarget, componentNoTarget);
}


}   // namespace

#include "mesh/mapping_between_meshes/manager/02_manager_initialize_composite.h"

namespace MappingBetweenMeshes
{

template<typename FunctionSpace1Type, typename FunctionSpace2Type>
void ManagerInitializeComposite::
initializeMappingsBetweenMeshes(const std::shared_ptr<FunctionSpace1Type> functionSpace1,
                                const std::shared_ptr<FunctionSpace2Type> functionSpace2)
{
  InitializeMappingsHelper<FunctionSpace1Type,FunctionSpace2Type>::
    initializeMappingsBetweenMeshes(shared_from_this(), functionSpace1, functionSpace2);
}

template<int D, typename BasisFunctionType, typename FunctionSpace2Type>
void InitializeMappingsHelper<
  FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,
  FunctionSpace2Type
>::initializeMappingsBetweenMeshes(std::shared_ptr<ManagerInitializeComposite> mappingsManager,
                                   const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>> functionSpace1,
                                   const std::shared_ptr<FunctionSpace2Type> functionSpace2)
{
  using SubFunctionSpaceType = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>;

  // add log message, to be included in the log file
  std::stringstream logMessage;
  logMessage << "  Mesh (\"" << functionSpace1->meshName() << "\") is composite and has "
    << functionSpace1->subFunctionSpaces().size() << " submeshes,";
  for (std::shared_ptr<SubFunctionSpaceType> &subFunctionSpace : functionSpace1->subFunctionSpaces())
    logMessage << " \"" << subFunctionSpace->meshName() << "\"";
  logMessage << ".";

  mappingsManager->addLogMessage(logMessage.str());

  // initialize mapping from all sub meshes

  // loop over sub function spaces
  for (std::shared_ptr<SubFunctionSpaceType> &subFunctionSpace : functionSpace1->subFunctionSpaces())
  {
    InitializeMappingsHelper<SubFunctionSpaceType,FunctionSpace2Type>::
      initializeMappingsBetweenMeshes(subFunctionSpace, functionSpace2);
  }

  mappingsManager->initializeMappingsBetweenMeshesFromSettings(functionSpace1, functionSpace2);
}

template<typename FunctionSpace1Type, typename FunctionSpace2Type>
void InitializeMappingsHelper<FunctionSpace1Type,FunctionSpace2Type>::
initializeMappingsBetweenMeshes(std::shared_ptr<ManagerInitializeComposite> mappingsManager,
                                const std::shared_ptr<FunctionSpace1Type> functionSpace1,
                                const std::shared_ptr<FunctionSpace2Type> functionSpace2)
{
  LOG(DEBUG) << "initializeMappingsBetweenMeshes source: \"" << functionSpace1->meshName() << "\" target: \"" << functionSpace2->meshName() << "\".";

  // add log message, to be included in the log file
  std::stringstream logMessage;
  logMessage << "  Initialize mappings between meshes \"" << functionSpace1->meshName()
    << "\", type " << StringUtility::demangle(typeid(FunctionSpace1Type).name())
    << "\n                                 and \"" << functionSpace2->meshName()
    << "\", type " << StringUtility::demangle(typeid(FunctionSpace2Type).name());

  mappingsManager->addLogMessage(logMessage.str());

  // this method is called in initialize of the operator splitting and should create the mapping objects between functionSpace1 and functionSpace2
  assert(FunctionSpace1Type::dim() <= FunctionSpace2Type::dim());   // the first mesh has to be lower dimensional than the second (or equal).

  mappingsManager->initializeMappingsBetweenMeshesFromSettings(functionSpace1, functionSpace2);
}

}   // namespace

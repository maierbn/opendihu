#pragma once

#include <Python.h>  // has to be the first included header
#include <map>
#include <vector>

#include "mesh/mapping_between_meshes/mapping/01_implementation.h"
#include "mesh/mapping_between_meshes/mapping/02_composite.h"
#include "field_variable/00_field_variable_base.h"

namespace MappingBetweenMeshes
{
  
class ManagerImplementation
{
public:
  //! constructor
  ManagerImplementation(PythonConfig specificSettings);

  //! get the mapping from source mesh to target mesh, create it if it does not yet exist
  template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  std::shared_ptr<MappingBetweenMeshes<typename FunctionSpaceSourceType::FunctionSpace, typename FunctionSpaceTargetType::FunctionSpace>>
  mappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                       std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget);

  //! prepare a mapping to the fieldVariableTarget, this zeros the targetFactorsSum for the field variable
  template<typename FieldVariableTargetType>
  void prepareMappingLowToHigh(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget);

  //! map data from on component of the source field variable to one component of the target field variable.
  //! Dimensionality source function space <= dimensionality target function space (e.g. 1D -> 3D)
  //! This has to be called between prepareMapping and finalizeMapping, can be called multiple times with different source meshes.
  //! If componentNoSource and componentNoTarget are both -1, map the whole field variable with all components
  template<typename FieldVariableSourceType, typename FieldVariableTargetType>
  void mapLowToHighDimension(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
                             std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget = -1);

  //! map complete data (all components) from source field variable to the target field variable.
  //! Dimensionality source function space >= dimensionality target function space (e.g. 3D -> 1D)
  //! If componentNoSource and componentNoTarget are both -1, map the whole field variable with all components
  template<typename FieldVariableSourceType, typename FieldVariableTargetType>
  void mapHighToLowDimension(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
                             std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget = -1);

  //! finalize the mapping to the fieldVariableTarget, for a mapping of a single component, this computes the final values at the dofs from the accumulated values by dividing by the targetFactorSums
  //! if componentNoTarget == -1, this means all components were transferred
  template<typename FieldVariableTargetType>
  void finalizeMappingLowToHigh(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget);

  //! finalize the mapping to the fieldVariableTarget, for a mapping of all components, this computes the final values at the dofs from the accumulated values by dividing by the targetFactorSums
  template<typename FieldVariableTargetType>
  void finalizeMappingLowToHigh(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget);

protected:

  //! create MappingBetweenMeshes objects from the config and store them under mappingsBetweenMeshes_
  void storeMappingsBetweenMeshes();

  //! create a single MappingBetweenMeshes object, targetMeshPy still needs to be parsed
  void storeMappingBetweenMeshes(std::string sourceMeshName, PyObject *targetMeshPy);

  //! indicate that on the mesh with name, "initialize()" has been called and now check if mappingsBetweenMeshes_ can be initialized
  void checkInitializeMappingBetweenMeshes(std::string name);

  //! Add and initialize a mapping between meshes, this can be called in the code in initialize, when it is clear, that this mapping will be needed
  //! If it is not clear, whether it will be needed, call initializeMappingsBetweenMeshes instead.
  template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  std::shared_ptr<MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>>
    createMappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget);

  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key, for meshManager

  struct MappingWithSettings
  {
    std::shared_ptr<MappingBetweenMeshesBase> mapping;   //< the actual mapping between the two meshes
    double xiTolerance;         //< the xiTolerance setting, 0 for disabled
    bool enableWarnings;        //< if warnings should be shown if source dofs are outside the target mesh with the given xi tolerance
    bool compositeUseOnlyInitializedMappings;   //< if for composite source meshes the mapping should be created from also defined mappings from the sub meshes
  };

  std::map<std::string, std::shared_ptr<FieldVariable::FieldVariableBase>> targetFactorSum_;    ///< for every target function space the field variable which accumulates the interpolation factors
  std::map<std::string, std::map<std::string, MappingWithSettings>> mappingsBetweenMeshes_;  ///<["key mesh from"]["key mesh to"] mapping between meshes
};

}  // namespace

#include "mesh/mapping_between_meshes/manager/00_manager_implementation.tpp"

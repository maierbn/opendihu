#pragma once

#include <Python.h>  // has to be the first included header
#include <map>
#include <vector>

#include "mesh/mapping_between_meshes.h"
#include "field_variable/00_field_variable_base.h"

namespace Mesh
{
/**
 * This class creates and stores all used meshes.
 * Each mesh can be defined in the python config under "Meshes" with a name and other properties.
 * Various components of the program can later
 * request their mesh by a call to mesh(name).
 * If a mesh was not defined earlier, it is created on the fly when it is requested.
 */
class MappingBetweenMeshesManager
{
public:
  //! constructor
  MappingBetweenMeshesManager(PythonConfig specificSettings);

  //! check if a MappingBetweenMeshes object need to be created and initialized
  template<typename FunctionSpace1Type, typename FunctionSpace2Type>
  void initializeMappingsBetweenMeshes(const std::shared_ptr<FunctionSpace1Type> functionSpace1, const std::shared_ptr<FunctionSpace2Type> functionSpace2);

  //! prepare a mapping to the fieldVariableTarget, this zeros the targetFactorsSum for the field variable
  template<typename FieldVariableTargetType>
  void prepareMapping(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget);

  //! map data from the source to the target field variable. This has to be called between prepareMapping and finalizeMapping, can be called multiple times with different source meshes.
  template<typename FieldVariableSourceType, typename FieldVariableTargetType>
  void map(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
           std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget);

  //! finalize the mapping to the fieldVariableTarget, this computes the final values at the dofs from the accumulated values by dividing by the targetFactorSums
  template<typename FieldVariableTargetType>
  void finalizeMapping(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget);

protected:

  //! create MappingBetweenMeshes objects from the config and store them under mappingsBetweenMeshes_
  void storeMappingsBetweenMeshes();

  //! indicate that on the mesh with name, "initialize()" has been called and now check if mappingsBetweenMeshes_ can be initialized
  void checkInitializeMappingBetweenMeshes(std::string name);

  //! get the mapping from source mesh to target mesh
  std::shared_ptr<MappingBetweenMeshesBase> mappingBetweenMeshes(std::string sourceMeshName, std::string targetMeshName);

  //! add and initialize a mapping between meshes
  template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  std::shared_ptr<MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>>
    createMappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget);

  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key, for meshManager

  std::map<std::string, std::shared_ptr<FieldVariable::FieldVariableBase>> targetFactorSum_;    ///< for every target function space the field variable which accumulates the interpolation factors
  std::map<std::string, std::map<std::string, std::shared_ptr<MappingBetweenMeshesBase>>> mappingsBetweenMeshes_;  ///<["key mesh from"]["key mesh to"] mapping between meshes
};

}  // namespace

#include "mesh/mapping_between_meshes_manager.tpp"

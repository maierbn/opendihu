#pragma once

#include <Python.h>  // has to be the first included header
#include <map>
#include <vector>

#include "function_space/function_space.h"

#include "mesh/mapping_between_meshes/manager/00_manager_log.h"
#include "mesh/mapping_between_meshes/mapping/02_composite.h"
#include "field_variable/00_field_variable_base.h"

namespace MappingBetweenMeshes
{

class ManagerInitialize : public ManagerLog
{
public:
  //! constructor
  ManagerInitialize(PythonConfig specificSettings);

  //! get the mapping from source mesh to target mesh, create it if it does not yet exist
  template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  std::shared_ptr<MappingBetweenMeshes<typename FunctionSpaceSourceType::FunctionSpace, typename FunctionSpaceTargetType::FunctionSpace>>
  mappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                       std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget);

  //! check if the mapping from source to target mesh exists
  template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  bool hasMappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource,
                               std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget);

  //! helper function that call the initialization depending on which mapping direction has been set in the python config
  template<typename FunctionSpace1Type, typename FunctionSpace2Type>
  void initializeMappingsBetweenMeshesFromSettings(const std::shared_ptr<FunctionSpace1Type> functionSpace1,
                                                   const std::shared_ptr<FunctionSpace2Type> functionSpace2);

protected:

  //! create MappingBetweenMeshes objects from the config and store them under mappingsBetweenMeshes_
  void storeMappingsBetweenMeshes(PythonConfig specificSettings);

  //! create a single MappingBetweenMeshes object, targetMeshPy still needs to be parsed
  void storeMappingBetweenMeshes(std::string sourceMeshName, PyObject *targetMeshPy);

  //! Add and initialize a mapping between meshes, this can be called in the code in initialize, when it is clear, that this mapping will be needed
  //! If it is not clear, whether it will be needed, call initializeMappingsBetweenMeshes instead.
  template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  std::shared_ptr<MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>>
    createMappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget);

  PythonConfig specificSettings_;    //< python object containing the value of the python config dict with corresponding key, for meshManager

  struct MappingWithSettings
  {
    std::shared_ptr<MappingBetweenMeshesBase> mapping;   //< the actual mapping between the two meshes
    double xiTolerance;         //< the xiTolerance setting, 0 for disabled
    bool enableWarnings;        //< if warnings should be shown if source dofs are outside the target mesh with the given xi tolerance
    bool compositeUseOnlyInitializedMappings;   //< if for composite source meshes the mapping should be created from also defined mappings from the sub meshes
    bool isEnabledFixUnmappedDofs;              //< if the unmapped dofs in the target mesh should be fixed by interpolating in the source mesh
    double defaultValue;        //< default value that is used if a target dof has no source dof that provides any values
  };

  std::map<std::string, std::shared_ptr<FieldVariable::FieldVariableBase>> targetFactorSum_;  //< for every target function space the field variable which accumulates the interpolation factors
  std::map<std::string, double> defaultValues_;    //< for every target mesh name the default value to set all target field variables to where there are no source dofs

  std::map<std::string, std::map<std::string, MappingWithSettings>> mappingsBetweenMeshes_;   //<["key mesh from"]["key mesh to"] mapping between meshes

};

}  // namespace

#include "mesh/mapping_between_meshes/manager/01_manager_initialize.tpp"

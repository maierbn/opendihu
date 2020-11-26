#pragma once

#include <Python.h>  // has to be the first included header

#include <memory>
#include "control/types.h"

#include "mesh/mapping_between_meshes/mapping/00_construct.h"

namespace MappingBetweenMeshes
{

/**
 * This is a mapping between two meshes, e.g. one 1D fiber mesh and one 3D mesh.
 * The mapping mapLowToHighDimension is from source mesh (lower dimensionality) to target mesh (higher dimensionality).
 * The mapping mapHighToLowDimension is from FunctionSpaceTargetType to FunctionSpaceSourceType.
 * Also read the more detailed description of the MappingBetweenMeshes::Manager class.
 */
template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
class MappingBetweenMeshesImplementation :
  public MappingBetweenMeshesConstruct<FunctionSpaceSourceType,FunctionSpaceTargetType>
{
public:

  //! constructor
  using MappingBetweenMeshesConstruct<FunctionSpaceSourceType,FunctionSpaceTargetType>::MappingBetweenMeshesConstruct;

  //! map data between a single component of the field variables in the source and target function spaces
  template<int nComponentsSource, int nComponentsTarget>
  void mapLowToHighDimension(FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponentsSource> &fieldVariableSource, int componentNoSource,
                             FieldVariable::FieldVariable<FunctionSpaceTargetType,nComponentsTarget> &fieldVariableTarget, int componentNoTarget,
                             FieldVariable::FieldVariable<FunctionSpaceTargetType,1> &targetFactorSum);

  //! map data between all components of the field variables in the source and target function spaces
  template<int nComponents>
  void mapLowToHighDimension(FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponents> &fieldVariableSource,
                             FieldVariable::FieldVariable<FunctionSpaceTargetType,nComponents> &fieldVariableTarget,
                             FieldVariable::FieldVariable<FunctionSpaceTargetType,1> &targetFactorSum);

  //! map data between all components of the field variables in the source and target function spaces
  //! the naming of the parameters is such that source->target, however the data types have swapped source/target semantics, because the mapping is the reverse one.
  template<int nComponents>
  void mapHighToLowDimension(FieldVariable::FieldVariable<FunctionSpaceTargetType,nComponents> &fieldVariableSource,
                             FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponents> &fieldVariableTarget);

  //! map data between between a single component of the field variables in the source and target function spaces
  //! the naming of the parameters is such that source->target, however the data types have swapped source/target semantics, because the mapping is the reverse one.
  template<int nComponentsTarget, int nComponentsSource>
  void mapHighToLowDimension(FieldVariable::FieldVariable<FunctionSpaceTargetType,nComponentsTarget> &fieldVariableSource, int componentNoSource,
                             FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponentsSource> &fieldVariableTarget, int componentNoTarget);
};

}  // namespace

// the include is contained in mesh/mapping_between_meshes/manager/04_manager.h
//#include "mesh/mapping_between_meshes/mapping/01_implementation.tpp"

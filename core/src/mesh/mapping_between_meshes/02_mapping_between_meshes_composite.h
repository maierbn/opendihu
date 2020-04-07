#pragma once

#include <Python.h>  // has to be the first included header

#include "mesh/mapping_between_meshes/01_mapping_between_meshes_implementation.h"
#include "mesh/composite.h"

namespace Mesh
{

/** Normal class for any mesh, except composite meshes.
 */
template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
class MappingBetweenMeshes : 
  public MappingBetweenMeshesImplementation<FunctionSpaceSourceType,FunctionSpaceTargetType>
{
public:

  using MappingBetweenMeshesImplementation<FunctionSpaceSourceType,FunctionSpaceTargetType>::MappingBetweenMeshesImplementation;
};

/** Partial specialization for composite source mesh
 */
template<int D, typename BasisFunctionType, typename FunctionSpaceTargetType>
class MappingBetweenMeshes<FunctionSpace::FunctionSpace<CompositeOfDimension<D>,BasisFunctionType>, FunctionSpaceTargetType> : 
  public MappingBetweenMeshesImplementation<FunctionSpace::FunctionSpace<CompositeOfDimension<D>,BasisFunctionType>, FunctionSpaceTargetType>
{
public:

  typedef FunctionSpace::FunctionSpace<CompositeOfDimension<D>,BasisFunctionType> FunctionSpaceSourceType;

  //! constructor, the function spaces need to be initialized
  MappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget, 
                       double xiTolerance=0, bool enableWarnings=true, bool compositeUseOnlyInitializedMappings=false);
};

}  // namespace


#pragma once

#include <Python.h>  // has to be the first included header

#include "mesh/mapping_between_meshes/mapping/01_implementation.h"
#include "mesh/composite.h"

namespace MappingBetweenMeshes
{

/**
 * This is a mapping between two meshes, e.g. one 1D fiber mesh and one 3D mesh.
 * The mapping mapLowToHighDimension is from source mesh (lower dimensionality) to target mesh (higher dimensionality).
 * The mapping mapHighToLowDimension is from FunctionSpaceTargetType to FunctionSpaceSourceType.
 * Also read the more detailed description of the MappingBetweenMeshes::Manager class.
 */
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
class MappingBetweenMeshes<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>, FunctionSpaceTargetType> : 
  public MappingBetweenMeshesImplementation<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>, FunctionSpaceTargetType>
{
public:

  typedef FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType> FunctionSpaceSourceType;

  //! constructor, the function spaces need to be initialized
  MappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget, 
                                double xiTolerance=0, bool enableWarnings=true, bool compositeUseOnlyInitializedMappings=false);
};

}  // namespace

// the following include is in mapping_between_meshes/manager/02_manager.h to prevent cyclic dependencies
//#include "mesh/mapping_between_meshes/mapping/02_composite.tpp"


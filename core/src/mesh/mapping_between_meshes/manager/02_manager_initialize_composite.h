#pragma once

#include <Python.h>  // has to be the first included header
#include <map>
#include <vector>

#include "mesh/mapping_between_meshes/manager/01_manager_initialize.h"

namespace MappingBetweenMeshes
{

/** Implementation of initialization for composite meshes
  */
class ManagerInitializeComposite :
  public ManagerInitialize,
  public std::enable_shared_from_this<ManagerInitializeComposite>
{
public:
  //! constructor
  using ManagerInitialize::ManagerInitialize;

  //! check if a MappingBetweenMeshes object need to be created and initialized, then initialize it if needed
  template<typename FunctionSpace1Type, typename FunctionSpace2Type>
  void initializeMappingsBetweenMeshes(const std::shared_ptr<FunctionSpace1Type> functionSpace1,
                                       const std::shared_ptr<FunctionSpace2Type> functionSpace2);
};

/** This implements the method initializeMappingsBetweenMeshes,
 *  It initializes the mappings between two function spaces FunctionSpace1Type and FunctionSpace2Type.
 *
 *  This needs to be a helper class, because if FunctionSpace1Type is a composite mesh, it first needs
 *  to initialize all mappings between the submeshes of the composite mesh and the target mesh.
 */
template<typename FunctionSpace1Type, typename FunctionSpace2Type>
class InitializeMappingsHelper
{
public:

  //! version for non-composite mesh in functionSpace1
  static void initializeMappingsBetweenMeshes(std::shared_ptr<ManagerInitializeComposite> mappingsManager,
                                              const std::shared_ptr<FunctionSpace1Type> functionSpace1,
                                              const std::shared_ptr<FunctionSpace2Type> functionSpace2);
};

/** partial specialization for a composite mesh
 */
template<int D, typename BasisFunctionType, typename FunctionSpace2Type>
class InitializeMappingsHelper<
  FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,
  FunctionSpace2Type
>
{
public:

  //! version for composite mesh in functionSpace1
  static void initializeMappingsBetweenMeshes(std::shared_ptr<ManagerInitializeComposite> mappingsManager,
                                              const std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>> functionSpace1,
                                              const std::shared_ptr<FunctionSpace2Type> functionSpace2);
};

}  // namespace

#include "mesh/mapping_between_meshes/manager/02_manager_initialize_composite.tpp"

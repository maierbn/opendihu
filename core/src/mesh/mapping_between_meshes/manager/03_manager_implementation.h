#pragma once

#include <Python.h>  // has to be the first included header
#include <map>
#include <vector>

#include "mesh/mapping_between_meshes/manager/02_manager_initialize_composite.h"

namespace MappingBetweenMeshes
{
  
class ManagerImplementation : public ManagerInitializeComposite
{
public:
  //! constructor
  using ManagerInitializeComposite::ManagerInitializeComposite;

  //! prepare a mapping to the fieldVariableTarget, this zeros the targetFactorsSum for the field variable
  template<typename FieldVariableTargetType>
  void prepareMappingLowToHigh(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget);

  //! map data from on component of the source field variable to one component of the target field variable.
  //! Dimensionality source function space <= dimensionality target function space (e.g. 1D -> 3D)
  //! This has to be called between prepareMapping and finalizeMapping, can be called multiple times with different source meshes.
  //! If componentNoSource and componentNoTarget are both -1, map the whole field variable with all components
  template<typename FieldVariableSourceType, typename FieldVariableTargetType>
  void mapLowToHighDimension(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
                             std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget = -1);

  //! map complete data (all components) from source field variable to the target field variable.
  //! Dimensionality source function space >= dimensionality target function space (e.g. 3D -> 1D)
  //! Internally, this uses the mapping from target to source function space.
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

  //! repair invalid geometry values for 1D fibers after mapping from a 3D mesh
  template<typename FieldVariableTargetType>
  void repairMappedGeometryFibers(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget);

protected:

  //! set the component of the target field variable to all zero (if it is a component), or all components (if componentNoTarget == -1)
  template<typename FieldVariableTargetType>
  void zeroTargetFieldVariable(std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget);
};

}  // namespace

#include "mesh/mapping_between_meshes/manager/03_manager_implementation.tpp"

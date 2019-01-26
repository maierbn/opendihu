#pragma once

#include <Python.h>  // has to be the first included header

namespace Mesh
{

class MappingBetweenMeshesBase{};

/**
 * This is a mapping between two meshes, e.g. one 1D fiber mesh and one 3D mesh. The mapping is from source mesh to target mesh.
 */
template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
class MappingBetweenMeshes : public MappingBetweenMeshesBase
{
public:

  //! constructor, the function spaces need to be initialized
  MappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget);

  //! map data between field variables in the source and target function spaces
  template<int nComponentsSource, int nComponentsTarget>
  void map(FieldVariable::FieldVariable<FunctionSpaceSourceType,nComponentsSource> &fieldVariableSource, int componentNoSource,
           FieldVariable::FieldVariable<FunctionSpaceTargetType,nComponentsTarget> &fieldVariableTarget, int componentNoTarget);

private:

  std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource_;   ///< the function space of the mesh from which to map data
  std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget_;   ///< the function space of the mesh to which to map data

  struct targetDof_t
  {
    element_no_t elementNoLocal;   //< local element no of the target element
    std::array<double,FunctionSpaceTargetType::nDofsPerElement()> scalingFactors;          //< factors for the dofs of the element with which to scale the value
  };

  std::vector<targetDof_t> targetMappingInfo_;  ///< [localDofNo source functionSpace] information where in the target to store the value from local dof No of the source
};

}  // namespace

#include "mesh/mapping_between_meshes.tpp"

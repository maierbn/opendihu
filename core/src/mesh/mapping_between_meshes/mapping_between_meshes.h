#pragma once

#include <Python.h>  // has to be the first included header

#include <memory>
#include "control/types.h"

//#include "field_variable/field_variable.h"
namespace FieldVariable
{
  template<typename FunctionSpace,int nComponents>
  class FieldVariable;
}

namespace Mesh
{

class MappingBetweenMeshesBase{};

/**
 * This is a mapping between two meshes, e.g. one 1D fiber mesh and one 3D mesh.
 * The mapping mapLowToHighDimension is from source mesh (lower dimensionality) to target mesh (higher dimensionality).
 * The mapping mapHighToLowDimension is from FunctionSpaceTargetType to FunctionSpaceSourceType.
 * Also read the more detailed description of the MappingBetweenMeshesManager class.
 */
template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
class MappingBetweenMeshes : public MappingBetweenMeshesBase
{
public:

  //! constructor, the function spaces need to be initialized
  MappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget, double xiTolerance=0);

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

private:

  std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource_;   ///< the function space of the mesh from which to map data
  std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget_;   ///< the function space of the mesh to which to map data

  struct targetDof_t
  {
    element_no_t elementNoLocal;   //< local element no of the target element (high dim)
    std::array<double,FunctionSpaceTargetType::nDofsPerElement()> scalingFactors;          //< factors for the dofs of the element with which to scale the value
    bool mapThisDof;    //< if this source dof should be mapped to the target dofs in elementNoLocal, if this is false, the dof is outside of the target mesh
  };

  std::vector<targetDof_t> targetMappingInfo_;  ///< [localDofNo source functionSpace (low dim)] information where in the target (high dim) to store the value from local dof No of the source (low dim)
  double maxAllowedXiTolerance_;    //< the maximum tolerance in the parameter space (threshold for xi) by which a node is considered still inside the element
};

}  // namespace

#include "mesh/mapping_between_meshes/mapping_between_meshes.tpp"

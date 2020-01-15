#pragma once

#include <Python.h>  // has to be the first included header
#include <map>
#include <vector>

#include "mesh/mapping_between_meshes/mapping_between_meshes.h"
#include "field_variable/00_field_variable_base.h"

namespace Mesh
{
/**
 * This class manages data mappings between different meshes, possibly of different dimensionalities, e.g. 1D->3D and inverse.
 * One mapping object handles both directions of the mapping, e.g. 1D->3D and 3D->1D. But the method to call is different, see below.
 *
 * At first, initializeMappingsBetweenMeshes has to be called once. (This is done at initalization of operator splitting.)
 * Here, the lower dimensional mesh has to be the first mesh and the higher dimensional mesh the second mesh.
 * Then, to perform the mapping, start with prepareMapping, then map and then finalizeMapping.
 *
 * A linear mapping of data from source mesh/function space to target mesh/function space is performed.
 * There are the mentioned two different cases, one for mapping from a lower dimension mesh to a higher dimension mesh
 * and one for the opposite direction.
 * The high to low dimension case simply evaluates the low dimension mesh dofs at the locations of the high dimension dofs.
 * This is performed using the normal interpolation defined by the function space, e.g. linear for linear Lagrange elements in the higher dimensional function space.
 * The mapping from a low dimension mesh to a high dimension mesh works as follows:
 * For every dof in the source function space (lower dimension) their contribution to dofs in the target function space (higher dimension) are determined.
 * This is done using the elements of the target function space (higher dimension), the contribution factors are just the linear interpolation factors,
 * i.e. the evaluated basis functions of the higher dimensional function space.
 * The source values multiplied with these contribution factors are summed up at the target dofs. At the end, the values are divided by the
 * summed up contribution factors at the target dofs.
 *
 * There are also two versions, one where a single component out of the source field variable is transferred to the target field variable
 * (e.g. transfer only one state like Vm from Monodomain eq. to EMG solver)
 * and the one where all components are transferred (e.g. 3D geometry field from elasticity solver to fibers).
 *
 *
 * Example in 1D for mapLowToHighDimension (in this case source and target have both the same :
 *
 * meshes:      t0--s0--------t1-----s1-----t2       t=target dof position, s=source dot position, a source mesh with 2 dofs in mapped to a target mesh with 3 dofs
 * src values:       2                 4
 *
 * determine targetFactors (factors of contribution):
 *              0.8*s0------ 0.2*s0 -------- 0.5*s1
 *                          +0.5*s1
 * total:       0.8 --------- 0.7 ----------- 0.5
 *
 * accumulate values:
 *              1.6 --------- 2.4 ----------- 2.0
 *
 * divide by targetFactors:
 *              1.6/0.8 ------ 2.4/0.7 ------ 2.0/0.5
 * result:       2.0 ---------  3.43 --------  4.0
 *
 *
 *
 * Example for mapLowToHighDimension:
 * meshes:      s0--t0--------s1-----t1-----s2       t=target dof position, s=source dot position, a source mesh with 2 dofs in mapped to a target mesh with 3 dofs
 * src values:  2             3.43          4
 *
 * target:        +0.8*s0           0.5*s1
 *                +0.2*s1          +0.5*s2
 * result:        2.286              3.715
 *
 *
 * How to call the mapping methods:
 * Call initializeMappingsBetweenMeshes once for initialization.
 * 1. mapping from lower to higher dimension function space (e.g. 1D->3D):
 *   fieldVariableSource = lower dimension (e.g. 1D)
 *   fieldVariableTarget = higher dimension (e.g. 3D)
 *
 *   prepareMappingLowToHigh(fieldVariableTarget)
 *
 *   mapLowToHighDimension(fieldVariableSource, fieldVariableTarget)   or   mapLowToHighDimension(fieldVariableSource, componentNoSource, fieldVariableTarget, componentNoTarget)
 *   (can be repeated for multiple lower dimension meshes, that all together should be mapped to the higher dimension space, e.g. multiple 1D fibers to one 3D mesh)
 *
 *   finalizeMappingLowToHigh(fieldVariableTarget)   or   finalizeMapping(fieldVariableTarget, componentNoTarget)
 *
 * 2. mapping from higher to lower dimension fucntion space (e.g. 3D->1D)
 *   fieldVariableSource = higher dimension (e.g. 3D)
 *   fieldVariableTarget = lower dimension (e.g. 1D)
 *
 *   mapHighToLowDimension(fieldVariableSource, fieldVariableTarget)
 *
 * Note, how prepareMapping and finalizeMapping are only needed for case 1. Ghosts are handled correctly within prepareMapping and finalizeMapping.
 * The 2. case is completely local and does not set ghost dofs of fieldVariableTarget. If ghost values are needed afterwards, call fieldVariableTarget->startGhostManipulation().
 *
 * There is now also the possiblity to use the simplified methods "prepareMapping", "map" and "finalizeMapping".
 * No initialization is required. Just always call:
 *
 * prepareMapping(fieldVariableSource, fieldVariableTarget)
 * map(fieldVariableSource, componentNoSource, fieldVariableTarget, componentNoTarget)  // set both componentNos to -1 to map all components
 * finalizeMapping(fieldVariableSource, fieldVariableTarget)
 *
 *
 */
class MappingBetweenMeshesManager
{
public:
  //! constructor
  MappingBetweenMeshesManager(PythonConfig specificSettings);

  //! check if a MappingBetweenMeshes object need to be created and initialized
  //! FunctionSpace1Type should be the lower dimension function space, FunctionSpace2Type the higher dimension function space
  template<typename FunctionSpace1Type, typename FunctionSpace2Type>
  void initializeMappingsBetweenMeshes(const std::shared_ptr<FunctionSpace1Type> functionSpace1, const std::shared_ptr<FunctionSpace2Type> functionSpace2);

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

  //! Add and initialize a mapping between meshes, this can be called in the code in initialize, when it is clear, that this mapping will be needed
  //! If it is not clear, whether it will be needed, call initializeMappingsBetweenMeshes instead.
  template<typename FunctionSpaceSourceType, typename FunctionSpaceTargetType>
  std::shared_ptr<MappingBetweenMeshes<FunctionSpaceSourceType, FunctionSpaceTargetType>>
    createMappingBetweenMeshes(std::shared_ptr<FunctionSpaceSourceType> functionSpaceSource, std::shared_ptr<FunctionSpaceTargetType> functionSpaceTarget);

  //! Simplified methods that call the other methods

  //! prepare the mapping for meshes of any dimensionality, this can be called even if not needed
  template<typename FieldVariableSourceType, typename FieldVariableTargetType>
  void prepareMapping(std::shared_ptr<FieldVariableSourceType> fieldVariableSource,
                      std::shared_ptr<FieldVariableTargetType> fieldVariableTarget);

  //! map between meshes of any dimensionality,
  //! if the field variables share the same function space, do no mapping at all, but copy the values
  //! If componentNoSource and componentNoTarget are both -1, map the whole field variable with all components
  template<typename FieldVariableSourceType, typename FieldVariableTargetType>
  void map(std::shared_ptr<FieldVariableSourceType> fieldVariableSource, int componentNoSource,
           std::shared_ptr<FieldVariableTargetType> &fieldVariableTarget, int componentNoTarget = -1, bool avoidCopyIfPossible=true);

  //! finalize the mapping for meshes of any dimensionality, this can be called even if not needed
  //! If componentNoTarget is -1, this means all components were transferred
  template<typename FieldVariableSourceType, typename FieldVariableTargetType>
  void finalizeMapping(std::shared_ptr<FieldVariableSourceType> fieldVariableSource,
                       std::shared_ptr<FieldVariableTargetType> fieldVariableTarget, int componentNoTarget = -1);

protected:

  //! create MappingBetweenMeshes objects from the config and store them under mappingsBetweenMeshes_
  void storeMappingsBetweenMeshes();

  //! indicate that on the mesh with name, "initialize()" has been called and now check if mappingsBetweenMeshes_ can be initialized
  void checkInitializeMappingBetweenMeshes(std::string name);

  //! get the mapping from source mesh to target mesh
  std::shared_ptr<MappingBetweenMeshesBase> mappingBetweenMeshes(std::string sourceMeshName, std::string targetMeshName);

  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key, for meshManager

  struct MappingWithSettings
  {
    std::shared_ptr<MappingBetweenMeshesBase> mapping;   //< the actual mapping between the two meshes
    double xiTolerance;         //< the xiTolerance setting, 0 for disabled
  };

  std::map<std::string, std::shared_ptr<FieldVariable::FieldVariableBase>> targetFactorSum_;    ///< for every target function space the field variable which accumulates the interpolation factors
  std::map<std::string, std::map<std::string, MappingWithSettings>> mappingsBetweenMeshes_;  ///<["key mesh from"]["key mesh to"] mapping between meshes
};

}  // namespace

#include "mesh/mapping_between_meshes/mapping_between_meshes_manager.tpp"
#include "mesh/mapping_between_meshes/mapping_between_meshes_manager_convenience_methods.tpp"

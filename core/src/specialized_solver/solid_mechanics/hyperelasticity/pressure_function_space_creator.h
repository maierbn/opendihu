#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>

#include "mesh/structured_deformable.h"
#include "function_space/function_space.h"
#include "mesh/mesh_manager/mesh_manager.h"

namespace SpatialDiscretization
{

/** Helper class that creates a first-order pressure function space from the 2nd order displacements function space by only using the nodes on the corners of the quadratic elements
 */
template<typename MeshType>
class PressureFunctionSpaceCreator
{
};

/** Implementation for 3D structured deformable mesh
*/
template<>
class PressureFunctionSpaceCreator<Mesh::StructuredDeformableOfDimension<3>>
{
public:

  typedef Mesh::StructuredDeformableOfDimension<3> MeshType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunction::LagrangeOfOrder<2>> DisplacementsFunctionSpace;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunction::LagrangeOfOrder<1>> PressureFunctionSpace;

  //! create the first order discretized pressure function space from the second order displacements function space
  static std::shared_ptr<PressureFunctionSpace> createPressureFunctionSpace(std::shared_ptr<Mesh::Manager> meshManager, std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace);

  //! from a vector of values that correspond to a field variable on the displacements function space, extract the values that correspond to the pressure function space
  template<typename T>
  static void extractPressureFunctionSpaceValues(std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace, std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace,
                                                 const std::vector<T> &displacementsFunctionSpaceValues, std::vector<T> &pressureFunctionSpaceValues);
};

/** Implementation for 3D composite mesh
*/
template<>
class PressureFunctionSpaceCreator<Mesh::CompositeOfDimension<3>>
{
public:

  typedef Mesh::CompositeOfDimension<3> MeshType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunction::LagrangeOfOrder<2>> DisplacementsFunctionSpace;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunction::LagrangeOfOrder<1>> PressureFunctionSpace;

  // create the first order discretized pressure function space from the second order displacements function space
  static std::shared_ptr<PressureFunctionSpace> createPressureFunctionSpace(std::shared_ptr<Mesh::Manager> meshManager, std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace);

  //! from a vector of values that correspond to a field variable on the displacements function space, extract the values that correspond to the pressure function space
  template<typename T>
  static void extractPressureFunctionSpaceValues(std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace, std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace,
                                                 const std::vector<T> &displacementsFunctionSpaceValues, std::vector<T> &pressureFunctionSpaceValues);
};

}  // namespace

#include "specialized_solver/solid_mechanics/hyperelasticity/pressure_function_space_creator.tpp"

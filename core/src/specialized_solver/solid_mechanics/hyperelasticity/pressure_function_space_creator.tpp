#include "specialized_solver/solid_mechanics/hyperelasticity/pressure_function_space_creator.h"

#include "control/types.h"
#include <vector>

namespace SpatialDiscretization
{

// structured 3D deformable mesh
template<typename T>
void PressureFunctionSpaceCreator<Mesh::StructuredDeformableOfDimension<3>>::
extractPressureFunctionSpaceValues(std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace, std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace,
                                   const std::vector<T> &displacementsFunctionSpaceValues, std::vector<T> &pressureFunctionSpaceValues)
{
  node_no_t nNodesLocal[3] = {
    displacementsFunctionSpace->meshPartition()->nNodesLocalWithGhosts(0),
    displacementsFunctionSpace->meshPartition()->nNodesLocalWithGhosts(1),
    displacementsFunctionSpace->meshPartition()->nNodesLocalWithGhosts(2)
  };

  int linearMeshIndex = pressureFunctionSpaceValues.size();   // append to previous values in vector
  pressureFunctionSpaceValues.resize(pressureFunctionSpaceValues.size() + pressureFunctionSpace->nNodesLocalWithGhosts());

  // loop over linear nodes in the quadratic mesh
  for (int k = 0; k < nNodesLocal[2]; k += 2)
  {
    for (int j = 0; j < nNodesLocal[1]; j += 2)
    {
      for (int i = 0; i < nNodesLocal[0]; i += 2, linearMeshIndex++)
      {
        int index = k*nNodesLocal[0]*nNodesLocal[1] + j*nNodesLocal[0] + i;

        pressureFunctionSpaceValues[linearMeshIndex] = displacementsFunctionSpaceValues[index];
      }
    }
  }
}

// composite 3D mesh
template<typename T>
void PressureFunctionSpaceCreator<Mesh::CompositeOfDimension<3>>::
extractPressureFunctionSpaceValues(std::shared_ptr<typename PressureFunctionSpaceCreator<Mesh::CompositeOfDimension<3>>::DisplacementsFunctionSpace> displacementsFunctionSpace,
                                   std::shared_ptr<typename PressureFunctionSpaceCreator<Mesh::CompositeOfDimension<3>>::PressureFunctionSpace> pressureFunctionSpace,
                                   const std::vector<T> &displacementsFunctionSpaceValues, std::vector<T> &pressureFunctionSpaceValues)
{
  // get the sub function spaces of the composite function spaces
  const std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>>> &displacementsSubFunctionSpaces
     = displacementsFunctionSpace->subFunctionSpaces();
  const std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>>> &pressureSubFunctionSpaces
     = pressureFunctionSpace->subFunctionSpaces();

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < displacementsSubFunctionSpaces.size(); subMeshNo++)
  {
    // get current sub meshes
    std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>> displacementsSubFunctionSpace = displacementsSubFunctionSpaces[subMeshNo];
    std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>> pressureSubFunctionSpace = pressureSubFunctionSpaces[subMeshNo];

    // call extractPressureFunctionSpaceValues on sub meshes
    PressureFunctionSpaceCreator<Mesh::StructuredDeformableOfDimension<3>>::extractPressureFunctionSpaceValues(displacementsSubFunctionSpace, pressureSubFunctionSpace, displacementsFunctionSpaceValues, pressureFunctionSpaceValues);
  }
}

}  // namespace

#include "specialized_solver/solid_mechanics/hyperelasticity/pressure_function_space_creator.h"

#include "control/types.h"
#include <vector>

namespace SpatialDiscretization
{

// structured 3D deformable mesh
std::shared_ptr<typename PressureFunctionSpaceCreator<Mesh::StructuredDeformableOfDimension<3>>::PressureFunctionSpace>
PressureFunctionSpaceCreator<Mesh::StructuredDeformableOfDimension<3>>::
createPressureFunctionSpace(std::shared_ptr<Mesh::Manager> meshManager, std::shared_ptr<typename PressureFunctionSpaceCreator<Mesh::StructuredDeformableOfDimension<3>>::DisplacementsFunctionSpace> displacementsFunctionSpace)
{
  // compose name of pressure function space
  std::stringstream name;
  name << displacementsFunctionSpace->meshName() << "_pressure";

  // if a function space with this name and type already exists, return the existing function space
  if (meshManager->hasFunctionSpaceOfType<PressureFunctionSpace>(name.str()))
  {
    return meshManager->functionSpace<PressureFunctionSpace>(name.str());
  }

  // create 3D function space with linear basis functions
  std::vector<Vec3> nodePositionsLinearMesh;

  // loop over nodes of quadratic function space and extract nodes for linear function space
  for (node_no_t nodeIndexZ = 0; nodeIndexZ < displacementsFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(2); nodeIndexZ++)
  {
    for (node_no_t nodeIndexY = 0; nodeIndexY < displacementsFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(1); nodeIndexY++)
    {
      for (node_no_t nodeIndexX = 0; nodeIndexX < displacementsFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(0); nodeIndexX++)
      {
        node_no_t nodeNoLocal = nodeIndexZ*displacementsFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(1)*displacementsFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(0)
          + nodeIndexY*displacementsFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(0) + nodeIndexX;

        dof_no_t dofNoLocal = nodeNoLocal;   // no Hermite, i.e. no multiple dofs per node

        if (nodeIndexX % 2 == 0 && nodeIndexY % 2 == 0 && nodeIndexZ % 2 == 0)
        {
          nodePositionsLinearMesh.push_back(displacementsFunctionSpace->geometryField().getValue(dofNoLocal));
        }
      }
    }
  }

  std::array<element_no_t,3> nElementsPerCoordinateDirection;
  std::array<int,3> nRanksPerCoordinateDirection;

  for (int i = 0; i < 3; i++)
  {
    nElementsPerCoordinateDirection[i] = displacementsFunctionSpace->meshPartition()->nElementsLocal(i);
    nRanksPerCoordinateDirection[i] = displacementsFunctionSpace->meshPartition()->nRanks(i);
  }

  // call mesh manager to create the new function space
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace;
  pressureFunctionSpace = meshManager->createFunctionSpace<PressureFunctionSpace>(
    name.str(), nodePositionsLinearMesh, nElementsPerCoordinateDirection, nRanksPerCoordinateDirection);

  return pressureFunctionSpace;
}

// composite 3D mesh
std::shared_ptr<typename PressureFunctionSpaceCreator<Mesh::CompositeOfDimension<3>>::PressureFunctionSpace>
PressureFunctionSpaceCreator<Mesh::CompositeOfDimension<3>>::
createPressureFunctionSpace(std::shared_ptr<Mesh::Manager> meshManager, std::shared_ptr<typename PressureFunctionSpaceCreator<Mesh::CompositeOfDimension<3>>::DisplacementsFunctionSpace> displacementsFunctionSpace)
{
  std::string name = "pressureFunctionSpace";

  // if a function space with this name and type already exists, return the existing function space
  if (meshManager->hasFunctionSpaceOfType<PressureFunctionSpace>(name))
  {
    return meshManager->functionSpace<PressureFunctionSpace>(name);
  }
  
  // get the sub function spaces of the composite displacements function space
  const std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>> >> &displacementsSubFunctionSpaces
     = displacementsFunctionSpace->subFunctionSpaces();

  // vector for pressure sub function spaces
  std::vector<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>>> pressureSubFunctionSpaces(displacementsSubFunctionSpaces.size());

  // loop over sub meshes and create pressure function spaces
  for (int subMeshNo = 0; subMeshNo < displacementsSubFunctionSpaces.size(); subMeshNo++)
  {
    pressureSubFunctionSpaces[subMeshNo] = PressureFunctionSpaceCreator<Mesh::StructuredDeformableOfDimension<3>>::createPressureFunctionSpace(meshManager, displacementsSubFunctionSpaces[subMeshNo]);
  }

  // create composite pressure function space
  std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace;
  pressureFunctionSpace = meshManager->createFunctionSpace<PressureFunctionSpace>(name, pressureSubFunctionSpaces);

  return pressureFunctionSpace;
}

}  // namespace

#pragma once

#include <Python.h>  // has to be the first included header

#include <memory>
#include "control/types.h"
#include "mesh/mapping_between_meshes/02_mapping_between_meshes_composite.h"

namespace Mesh
{

/**
 * This is a mapping between two meshes, e.g. one 1D fiber mesh and one 3D mesh.
 * The mapping mapLowToHighDimension is from source mesh (lower dimensionality) to target mesh (higher dimensionality).
 * The mapping mapHighToLowDimension is from FunctionSpaceTargetType to FunctionSpaceSourceType.
 * Also read the more detailed description of the MappingBetweenMeshesManager class.
 */

}  // namespace


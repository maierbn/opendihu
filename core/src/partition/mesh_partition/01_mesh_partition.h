#pragma once

#include <memory>
#include <petscdmda.h>

#include "partition/mesh_partition/00_mesh_partition_base.h"
#include "control/types.h"
#include "partition/rank_subset.h"
#include "mesh/type_traits.h"
#include "mesh/face_t.h"

namespace Partition
{
 
/** Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: starting with 0, up to total number including ghost elements
 */
template<typename FunctionSpaceType, typename DummyForTraits = typename FunctionSpaceType::Mesh>
class MeshPartition
{
};

}   // namespace

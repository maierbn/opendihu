#pragma once

namespace Partition
{

enum values_representation_t
{
  representationLocal,            // the local vector contains the valid data, the other vectors may be outdated
  representationGlobal,           // the global vector contains the valid data, the other vectors may be outdated
  representationContiguous,   // the contiguous vector contains the valid data, the other vectors may be outdated
  representationGlobalSub    // the global Vec is set as a subvector of the contiguos Vec, of the first component (only for 1 component vectors)
};

extern const char *valuesRepresentationString[16];

};  // namespace

#pragma once

namespace Partition
{

enum values_representation_t
{
  representationLocal,            // the local vector contains the valid data, the other vectors may be outdated
  representationGlobal,           // the global vector contains the valid data, the other vectors may be outdated
  representationContiguous,   // the contiguous vector contains the valid data, the other vectors may be outdated
  representationInvalid       // this field variable is not usable, because data has been extracted by extractComponentShared, call restoreExtractedComponent to make it usable again

};

extern const char *valuesRepresentationString[16];

};  // namespace

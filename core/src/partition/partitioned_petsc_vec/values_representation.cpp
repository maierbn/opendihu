#include "partition/partitioned_petsc_vec/values_representation.h"

namespace Partition
{

const char *valuesRepresentationString[16] =
{
  "local",
  "global",
  "contiguous",
  "invalid",
  "combined-local",
  "combined-global",
  "no-vector"
};

} // namespace

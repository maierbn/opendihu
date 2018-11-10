#include "basis_function/hermite.h"

namespace BasisFunction
{

constexpr int Hermite::nDofsPerNode()
{
  return 2;
}

constexpr int Hermite::nDofsPerBasis()
{
  return 4;
}

}; // namespace

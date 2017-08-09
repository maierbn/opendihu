#pragma once

#include "mesh/mesh.h"
#include "control/types.h"

namespace Mesh
{

template<int D>
class Regular : public Mesh
{
public:
  static element_idx_t nElements();
private:
};

}  // namespace

#include "mesh/regular.tpp"
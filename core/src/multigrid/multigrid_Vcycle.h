#pragma once

#include "multigrid/multigrid.h"

namespace Multigrid
{

template<typename FiniteElement1, typename FiniteElement2>
class multigrid_Vcycle :
  public Multigrid<FiniteElement1,FiniteElement2>
{
public:
  //! constructor
  multigrid_Vcycle(DihuContext context);

  
protected:
};

}  // namespace

#include "multigrid/multigrid_Vcycle.tpp"

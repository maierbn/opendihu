#include "multigrid/multigrid_Vcycle.h"

#include "utility/python_utility.h"
namespace Multigrid
{

template<typename FiniteElement1, typename FiniteElement2>
multigrid_Vcycle<FiniteElement1,FiniteElement2>::
multigrid_Vcycle(DihuContext context) :
  Multigrid<FiniteElement1,FiniteElement2>(context, "multigrid_Vcycle")
{
}
};    // namespace

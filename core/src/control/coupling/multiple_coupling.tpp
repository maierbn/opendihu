#include "control/coupling/multiple_coupling.h"

namespace Control
{

template<typename FirstNestedSolver, typename... MoreNestedSolvers>
MultipleCoupling<FirstNestedSolver, MoreNestedSolvers...>::
MultipleCoupling(DihuContext context) :
  NestedCouplingSchemes<FirstNestedSolver, MoreNestedSolvers...>::type
  (
    context["MultipleCoupling"],
    FirstNestedSolver(context["MultipleCoupling"]["Term1"]),      // construct single solver from settings under "Term1"
    NestedCouplingSchemes<MoreNestedSolvers...>::createScheme(context["MultipleCoupling"], 2)   // construct nested coupling schemes, starting from "Term2"
  )
{
}

}  // namespace

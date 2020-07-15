#include "control/coupling/nested_coupling_schemes.h"

namespace Control
{

template<typename FirstNestedSolver, typename ... MoreNestedSolvers>
typename NestedCouplingSchemes<FirstNestedSolver, MoreNestedSolvers...>::type NestedCouplingSchemes<FirstNestedSolver, MoreNestedSolvers...>::
createScheme(DihuContext context, int termNo)
{
  // create first solver for the Coupling scheme
  std::stringstream name;
  name << "Term" << termNo;
  //FirstNestedSolver &&firstNestedSolver = FirstNestedSolver(context[name.str()]);

  // recursively create second nested solver for the Coupling scheme
  //typename NestedCouplingSchemes<MoreNestedSolvers...>::type &&secondNestedSolver = NestedCouplingSchemes<MoreNestedSolvers...>::createScheme(context, termNo+1);

  return Coupling<FirstNestedSolver,typename NestedCouplingSchemes<MoreNestedSolvers...>::type>(context, FirstNestedSolver(context[name.str()]), NestedCouplingSchemes<MoreNestedSolvers...>::createScheme(context, termNo+1));
}

/** Partial specialization for recursion end.
 */
template<typename LastNestedSolver>
typename NestedCouplingSchemes<LastNestedSolver>::type NestedCouplingSchemes<LastNestedSolver>::
createScheme(DihuContext context, int termNo)
{
  std::stringstream name;
  name << "Term" << termNo;

  return LastNestedSolver(context[name.str()]);
}

}  // namespace

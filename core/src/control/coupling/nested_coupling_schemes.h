#pragma once

#include <Python.h>  // has to be the first included header

namespace Control
{

/** Helper struct that defines nested coupling schemes from a given list of solvers:
 *  NestedCouplingSchemes<A,B,C,D>::type ->
 *  Coupling<A, Coupling<B, Coupling<C,D>>
 */
template<typename... NestedSolvers>
struct NestedCouplingSchemes;

template<typename FirstNestedSolver, typename... MoreNestedSolvers>
struct NestedCouplingSchemes<FirstNestedSolver, MoreNestedSolvers...>
{
  typedef Coupling<FirstNestedSolver,typename NestedCouplingSchemes<MoreNestedSolvers...>::type> type;

  static type createScheme(DihuContext context, int termNo=1);
};

/** Partial specialization for recursion end.
 */
template<typename LastNestedSolver>
struct NestedCouplingSchemes<LastNestedSolver>
{
  typedef LastNestedSolver type;

  static type createScheme(DihuContext context, int termNo=1);
};

}  // namespace

#include "control/coupling/nested_coupling_schemes.tpp"

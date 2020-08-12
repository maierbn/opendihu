#pragma once

#include "control/coupling/coupling.h"
#include "control/coupling/nested_coupling_schemes.h"

namespace Control
{

/** MultipleCoupling provides a shortcut for multiple nested Coupling<A, Coupling<B, Coupling<C,D>>> schemes.
 *  Instead of nested settings, such as
 *  "Coupling": {
 *    "Term1": {...},
 *    "Term2": {
 *      "Coupling": {
 *        "Term1": {...},
 *        "Term2": {
 *          "Coupling": {
 *            "Term1": {...},
 *            "Term2": {...}
 *          }
 *        }
 *      }
 *    }
 *  },
 *  we have the terms on the same level, i.e.
 *  "MultipleCoupling": {
 *     "Term1": {...},
 *     "Term2": {...},
 *     "Term3": {...},
 *     "Term4": {...}
 *  }
 *
  */
template<typename FirstNestedSolver, typename... MoreNestedSolvers>
class MultipleCoupling :
  public NestedCouplingSchemes<FirstNestedSolver, MoreNestedSolvers...>::type    // inherit from outer-most Coupling<>
{
public:

  typedef typename NestedCouplingSchemes<FirstNestedSolver, MoreNestedSolvers...>::type NestedCouplingsType;  // type of the nested Coupling schemes

  //! constructor
  MultipleCoupling(DihuContext context);

};

}  // namespace

#include "control/coupling/multiple_coupling.tpp"

#pragma once

#include "equation/equations.h"
#include "basis_function/lagrange.h"
#include "mesh/structured_regular_fixed.h"

namespace Equation
{

// Equation of the form L = u_t
template<typename Term>
using usesTimeStepping = std::enable_if_t<Term::usesTimeStepping, Term>;

// Equations that include Δu
template<typename Term>
using hasLaplaceOperator = std::enable_if_t<Term::hasLaplaceOperator, Term>;

// Equations that include ∇•(A∇u)
template<typename Term>
using hasGeneralizedLaplaceOperator = std::enable_if_t<Term::hasGeneralizedLaplaceOperator, Term>;

// Equations that can have a non-zero rhs (Lu = f), but not Lu = u_t
template<typename Term>
using hasRhsNoTimestepping = std::enable_if_t<
  Term::hasRhs && !Term::usesTimeStepping,
  Term
>;

// Equations that can have a non-zero rhs (Lu = f)
template<typename Term>
using hasRhs = std::enable_if_t<Term::hasRhs,Term>;

// Equations that can have no rhs and not time stepping, i.e. Laplace (Lu = 0) or solid mechanics
template<typename Term>
using hasNoRhs = std::enable_if_t<!Term::hasRhs,Term>;

// Equations of solid mechanics
template<typename Term>
using isSolidMechanics = std::enable_if_t<Term::isSolidMechanics,Term>;

// compressible material
template<typename Term>
using isCompressible = std::enable_if_t<Term::isSolidMechanics && !Term::isIncompressible,Term>;

// incompressible material
template<typename Term>
using isIncompressible = std::enable_if_t<Term::isIncompressible,Term>;

template<typename BasisFunction, typename Mesh, typename Term>
using doesNotUseStencils = std::enable_if_t<
  // not linear Lagrange on regular fixed mesh
  !(std::is_same<BasisFunction, ::BasisFunction::LagrangeOfOrder<1>>::value  
    && (std::is_same<Mesh, ::Mesh::StructuredRegularFixedOfDimension<1>>::value
        || std::is_same<Mesh, ::Mesh::StructuredRegularFixedOfDimension<2>>::value
        || std::is_same<Mesh, ::Mesh::StructuredRegularFixedOfDimension<3>>::value))
  // or has a generalized laplace operator
  || Term::hasGeneralizedLaplaceOperator
  ,
  Mesh
>;

template<typename BasisFunction, typename Mesh, typename Term>
using doesNotUseStencilsNorSolidMechanics = std::enable_if_t<
  // not linear Lagrange on regular fixed mesh
  (!(std::is_same<BasisFunction, ::BasisFunction::LagrangeOfOrder<1>>::value  
    && (std::is_same<Mesh, ::Mesh::StructuredRegularFixedOfDimension<1>>::value
        || std::is_same<Mesh, ::Mesh::StructuredRegularFixedOfDimension<2>>::value
        || std::is_same<Mesh, ::Mesh::StructuredRegularFixedOfDimension<3>>::value))
  // or has a generalized laplace operator
  || Term::hasGeneralizedLaplaceOperator)
  
  // also not solid mechanics
  && !Term::isSolidMechanics
  ,
  Mesh
>;


} // namespace Equation

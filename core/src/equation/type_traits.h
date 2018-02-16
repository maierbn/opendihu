#pragma once

#include "equation/diffusion.h"
#include "equation/poisson.h"
#include "equation/laplace.h"
#include "equation/static.h"
#include "equation/solid_mechanics.h"
#include "equation/compressible_mooney_rivlin.h"

namespace Equation
{

// Equation of the form Δu = f
template<class Term>
using hasLaplaceOperatorWithRhs = std::enable_if_t<
  std::is_same<Term, Static::Poisson>::value
  || std::is_same<Term, Static::GeneralizedPoisson>::value,
  Term
>;

// Equation of the form Δu = u_t
template<class Term>
using hasLaplaceOperatorWithTimeStepping = std::enable_if_t<
  std::is_same<Term, Dynamic::Diffusion>::value,
  Term
>;

// Equations that include Δu
template<class Term>
using hasLaplaceOperator = std::enable_if_t<
  std::is_same<Term, Static::Laplace>::value
  || std::is_same<Term, Static::Poisson>::value
  || std::is_same<Term, Static::GeneralizedPoisson>::value
  || std::is_same<Term, Dynamic::Diffusion>::value,
  Term
>;

// Equations that include a proper rhs
template<class Term>
using hasRhs = std::enable_if_t<
  std::is_same<Term, Static::Poisson>::value
  || std::is_same<Term, Static::GeneralizedPoisson>::value
  || std::is_same<Term, Dynamic::Diffusion>::value,
  Term
>;

// Equations of solid mechanics
template<class Term>
using isSolidMechanics = std::enable_if_t<
  std::is_same<Term, Static::SolidMechanics>::value
  || std::is_same<Term, Static::CompressibleMooneyRivlin>::value,
  Term
>;



} // namespace Equation

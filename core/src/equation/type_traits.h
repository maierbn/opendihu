#pragma once

namespace Equation
{

// Equation of the form Δu = f
/*template<class Term>
using hasLaplaceOperatorWithRhs = std::enable_if_t<
  std::is_same<Term, Static::Poisson>
  || std::is_same<Term, Dynamic::Diffusion>,
  Term
>;
*/
template<class Term>
using hasLaplaceOperatorWithRhs = std::enable_if_t<
  true,
  Term
>;


} // namespace Equation
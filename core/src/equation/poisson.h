#pragma once

#include "equation/static.h"

namespace Equation
{
namespace Static
{

/** Δu = f
  */
class Poisson : public Static
{
};

/** ∇•A∇u = f
  */
class GeneralizedPoisson : public Static
{
};

}  // namespace
}  // namespace

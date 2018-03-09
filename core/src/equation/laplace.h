#pragma once

#include "equation/static.h"

namespace Equation
{
namespace Static
{

/** Δu = 0
  */
class Laplace : public Static
{
};

/** ∇•A∇u = 0
  */
class GeneralizedLaplace : public Static
{
};

}  // namespace
}  // namespace

#pragma once

#include "equation/dynamic.h"

namespace Equation
{
namespace Dynamic
{

/** Δu = f
  */
class IsotropicDiffusion : public Dynamic
{
};

/** ∇•A∇u = f
  */
class AnisotropicDiffusion : public Dynamic
{
};

}  // namespace
}  // namespace

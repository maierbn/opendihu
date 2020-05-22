#pragma once

#include <Python.h>  // has to be the first included header

#include <vector>
#include <Vc/Vc>

namespace SpatialDiscretization
{

template<typename T>
class ExpressionHelper {};

/** Helper class that inserts variables in a SEMT symbolic expression.
 * This partial specialization is for normal double values.
 */
template<>
class ExpressionHelper<double>
{
public:

  // apply the SEMT expression to the given variables
  template<typename SEMTExpressionType>
  static double apply(SEMTExpressionType &expression, const std::vector<double> &variables);
};

/** Partial specialization for Vc::double_v, i.e. vectorized apply for multiple sets of values at once
 */
template<>
class ExpressionHelper<Vc::double_v>
{
public:

  // apply the SEMT expression to the given variables
  template<typename SEMTExpressionType>
  static Vc::double_v apply(SEMTExpressionType &expression, const std::vector<Vc::double_v> &variables);
};

}  // namespace

#include "specialized_solver/solid_mechanics/hyperelasticity/expression_helper.tpp"

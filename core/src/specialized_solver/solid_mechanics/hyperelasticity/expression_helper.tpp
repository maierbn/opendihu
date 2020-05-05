#include "specialized_solver/solid_mechanics/hyperelasticity/expression_helper.h"

namespace SpatialDiscretization
{

template<typename SEMTExpressionType>
double ExpressionHelper<double>::apply(SEMTExpressionType &expression, const std::vector<double> &variables)
{
  return expression.apply(variables);
}

template<typename SEMTExpressionType>
Vc::double_v ExpressionHelper<Vc::double_v>::apply(SEMTExpressionType &expression, const std::vector<Vc::double_v> &variables)
{
  Vc::double_v result;
  std::vector<double> variablesVector(variables.size());

  // loop over the components of the vectorized data type
  for (int vcComponentNo = 0; vcComponentNo < Vc::double_v::size(); vcComponentNo++)
  {
    // loop over variables and set the variables vector for the current vc component
    for (int i = 0; i < variables.size(); i++)
    {
      variablesVector[i] = variables[i][vcComponentNo];
    }

    // apply expression for current component
    result[vcComponentNo] = expression.apply(variablesVector);
  }
  return result;
}

}  // namespace


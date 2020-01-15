#include "spatial_discretization/finite_element_method/integrand/integrand_stiffness_matrix_linear_elasticity.h"

#include <array>

#include "utility/math_utility.h"

namespace SpatialDiscretization
{

//integrand for stiffness matrix of finite elasticity
template<int D, typename EvaluationsType,typename FunctionSpaceType,typename Term>
EvaluationsType IntegrandStiffnessMatrix<D,EvaluationsType,FunctionSpaceType,D,Term,Equation::isLinearElasticity<Term>>::
evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,D,Term> &data, const std::array<Vec3,D> &jacobian,
                  element_no_t elementNoLocal, const std::array<double,D> xi)
{
  EvaluationsType evaluations{};

  double integrationFactor = MathUtility::computeIntegrationFactor<D>(jacobian);
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // compute gradient values in parameter space
  std::array<VecD<D>,nDofsPerElement> gradPhi = data.functionSpace()->getGradPhi(xi);

  // compute inverse jacobian
  std::array<Vec3,nDofsPerElement> geometryValues;
  data.functionSpace()->getElementGeometry(elementNoLocal, geometryValues);
  Tensor2<D> inverseJacobian = data.functionSpace()->getInverseJacobian(geometryValues, elementNoLocal, xi);


  // loop over entries (La, Mb) of stiffness matrix
  for (int indexL = 0; indexL < nDofsPerElement; indexL++)
  {
    for (int indexM = 0; indexM < nDofsPerElement; indexM++)
    {
      //! jacobianParameterSpace[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      //! inverseJacobianParameterSpace[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

      VecD<D> dPhiL_dX = inverseJacobian * gradPhi[indexL];
      VecD<D> dPhiM_dX = inverseJacobian * gradPhi[indexM];

      for (int a = 0; a < D; a++)
      {
        for (int b = 0; b < D; b++)
        {
          // loop over internal dimensions
          for (int c = 0; c < D; c++)
          {
            for (int d = 0; d < D; d++)
            {
              double integrand = data.linearStiffness(a, d, b, c) * dPhiL_dX[d] * dPhiM_dX[c] * integrationFactor;

              //if (a != b)
              //  LOG(DEBUG) << "abcd=" << a << b << c << d << ": stiffness=" << data.linearStiffness(a, d, b, c) << ", integrand=" << integrand
              //    << " at " << indexL << "," << indexM << " (" << indexL*D + a << "," << indexM*D + b << ")";

              if (!std::isfinite(integrand))
              {
                LOG(ERROR) << "Entry is not finite, abcd=" << a << b << c << d << ": stiffness=" << data.linearStiffness(a, d, b, c)
                  << ", integrand=" << integrand << " at " << indexL << "," << indexM << " (" << indexL*D + a << "," << indexM*D + b << ")";
              }

              VLOG(2) << "  dof pair (" << indexL << "," << indexM << "), "
                << "component (" << a << "," << b << ") internal (" << c << "," << d << "), integrated value: " << integrand;

              evaluations(indexL*D + a, indexM*D + b) += integrand;
            }
          }

          VLOG(2) << "  dof pair (" << indexL << "," << indexM << "), "
            << "component (" << a << "," << b << "), integrated value: " <<evaluations(indexL*D + a, indexM*D + b);
        }
      }
    }
  }

  return evaluations;
}

} // namespace

#include "specialized_solver/muscle_contraction_solver.h"

#include <omp.h>
#include <sstream>

#include "utility/math_utility.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

template<typename MeshType,typename Term,bool withLargeOutputFiles>
void MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
computeLambda()
{
  typedef typename DynamicHyperelasticitySolverType::HyperelasticitySolverType::DisplacementsFieldVariableType DisplacementsFieldVariableType;
  typedef typename DynamicHyperelasticitySolverType::HyperelasticitySolverType::Data::DeformationGradientFieldVariableType DeformationGradientFieldVariableType;
  typedef typename Data::ScalarFieldVariableType FieldVariableType;

  std::shared_ptr<DisplacementsFieldVariableType> fiberDirectionVariable;
  std::shared_ptr<DeformationGradientFieldVariableType> deformationGradientVariable;
  std::shared_ptr<DeformationGradientFieldVariableType> fDotVariable;
  std::shared_ptr<DisplacementsFieldVariableType> velocitiesVariable;

  if (isDynamic_)
  {
    fiberDirectionVariable = dynamicHyperelasticitySolver_->hyperelasticitySolver().data().fiberDirection();
    deformationGradientVariable = dynamicHyperelasticitySolver_->hyperelasticitySolver().data().deformationGradient();
    fDotVariable = dynamicHyperelasticitySolver_->hyperelasticitySolver().data().deformationGradientTimeDerivative();
    velocitiesVariable = dynamicHyperelasticitySolver_->data().velocities();
  }
  else
  {
    fiberDirectionVariable = staticHyperelasticitySolver_->data().fiberDirection();
    deformationGradientVariable = staticHyperelasticitySolver_->data().deformationGradient();
    fDotVariable = staticHyperelasticitySolver_->data().deformationGradientTimeDerivative();
    velocitiesVariable = staticHyperelasticitySolver_->data().velocities();
  }

  // compute lambda and \dot{lambda} (contraction velocity)
  //std::shared_ptr<DisplacementsFieldVariableType> displacementsVariable = dynamicHyperelasticitySolver_->data().displacements();

  std::shared_ptr<FieldVariableType> lambdaVariable = data_.lambda();
  std::shared_ptr<FieldVariableType> lambdaDotVariable = data_.lambdaDot();

  // loop over local degrees of freedom
  for (dof_no_t dofNoLocal = 0; dofNoLocal < data_.functionSpace()->nDofsLocalWithoutGhosts(); dofNoLocal++)
  {
    const Vec3 fiberDirection = fiberDirectionVariable->getValue(dofNoLocal);
    //const Vec3 displacement = displacementsVariable->getValue(dofNoLocal);
    const VecD<9> deformationGradientValues = deformationGradientVariable->getValue(dofNoLocal);
    const VecD<9> fDotValues = fDotVariable->getValue(dofNoLocal);

    // create matrix, deformationGradientValues are in row-major order
    MathUtility::Matrix<3,3> deformationGradient(deformationGradientValues);
    MathUtility::Matrix<3,3> fDot(fDotValues);

    // fiberDirection is normalized
    assert (MathUtility::norm<3>(fiberDirection) - 1.0 < 1e-10);

    // get deformation gradient, project lambda and lambda dot
    // dx = F dX, dx^2 = C dX^2
    // λ = ||dx•a0||
    // project displacements on normalized fiberDirection a0
    //const double lambda = displacement[0] * fiberDirection[0] + displacement[1] * fiberDirection[1] + displacement[2] * fiberDirection[2];

    // convert fiber direction from reference configuration into current configuration
    Vec3 fiberDirectionCurrentConfiguration = deformationGradient * fiberDirection;   // F a0
    // λ = ||F a0||
    const double lambda = MathUtility::norm<3>(fiberDirectionCurrentConfiguration);   // stretch in current configuration

    // exemplary derivative of λ for dim=2:
    //  λ = ||F a0|| = sqrt[(F11*a1 + F12*a2)^2 + (F21*a1 + F22*a2)^2]
    // d/dt ||F a0|| = 1/(2*sqrt[(F11*a1 + F12*a2)^2 + (F21*a1 + F22*a2)^2])
    //                 * (2*(F11*a1 + F12*a2)*(F11'*a1 + F12'*a2) + 2*(F21*a1 + F22*a2)*(F21'*a1 + F22'*a2))
    //               = (F a0) • (F' a0) / ||F a0||

    // compute lambda dot
    // d/dt λ = d/dt ||F a0|| = (F a0) • (Fdot a0) / ||F a0||   (where Fdot = d/dt F)
    // d/dt dx = d/dt F
    // d/dt lambda = d/dt ||dx•a0|| = 1 / ||Fa|| (F a0 • Fdot a0) = 1/lambda (Fa • Fdot a0)
    Vec3 FdotA0 = fDot * fiberDirection;
    const double lambdaDot = 1 / lambda * MathUtility::dot(fiberDirectionCurrentConfiguration, FdotA0) * lambdaDotScalingFactor_;

    lambdaVariable->setValue(dofNoLocal, lambda);
    lambdaDotVariable->setValue(dofNoLocal, lambdaDot);
  }

  lambdaVariable->zeroGhostBuffer();
  lambdaVariable->finishGhostManipulation();
  lambdaVariable->startGhostManipulation();

  lambdaDotVariable->zeroGhostBuffer();
  lambdaDotVariable->finishGhostManipulation();
  lambdaDotVariable->startGhostManipulation();
}

template<typename MeshType,typename Term,bool withLargeOutputFiles>
void MuscleContractionSolver<MeshType,Term,withLargeOutputFiles>::
computeActiveStress()
{
  LOG(DEBUG) << "computeActiveStress";

  typedef typename DynamicHyperelasticitySolverType::HyperelasticitySolverType::StressFieldVariableType StressFieldVariableType;
  typedef typename DynamicHyperelasticitySolverType::HyperelasticitySolverType::DisplacementsFieldVariableType DisplacementsFieldVariableType;
  typedef typename Data::ScalarFieldVariableType FieldVariableType;

  std::shared_ptr<StressFieldVariableType> activePK2StressVariable;
  std::shared_ptr<DisplacementsFieldVariableType> fiberDirectionVariable;

  if (isDynamic_)
  {
    activePK2StressVariable = dynamicHyperelasticitySolver_->hyperelasticitySolver().data().activePK2Stress();
    fiberDirectionVariable = dynamicHyperelasticitySolver_->hyperelasticitySolver().data().fiberDirection();
  }
  else
  {
    activePK2StressVariable = staticHyperelasticitySolver_->data().activePK2Stress();
    fiberDirectionVariable = staticHyperelasticitySolver_->data().fiberDirection();
  }

  std::shared_ptr<FieldVariableType> lambdaVariable = data_.lambda();
  std::shared_ptr<FieldVariableType> gammaVariable = data_.gamma();

  // Heidlauf 2013: "Modeling the Chemoelectromechanical Behavior of Skeletal Muscle Using the Parallel Open-Source Software Library OpenCMISS", p.4, Eq. (11)

  const double lambdaOpt = 1.2;

  // loop over local degrees of freedom
  for (dof_no_t dofNoLocal = 0; dofNoLocal < data_.functionSpace()->nDofsLocalWithoutGhosts(); dofNoLocal++)
  {
    const Vec3 fiberDirection = fiberDirectionVariable->getValue(dofNoLocal);

    const double lambda = lambdaVariable->getValue(dofNoLocal);
    const double gamma = gammaVariable->getValue(dofNoLocal);
    const double lambdaRelative = lambda / lambdaOpt;

    // compute f function
    double f = 1.0;

    if (enableForceLengthRelation_)
    {
      if (0.6 <= lambdaRelative && lambdaRelative <= 1.4)
      {
        f = -25./4 * lambdaRelative*lambdaRelative + 25./2 * lambdaRelative - 5.25;
      }
    }

    const double factor = 1./lambda * pmax_ * f * gamma;

    // Voigt notation:
    // [0][0] -> [0];
    // [1][1] -> [1];
    // [2][2] -> [2];
    // [0][1] -> [3];
    // [1][0] -> [3];
    // [1][2] -> [4];
    // [2][1] -> [4];
    // [0][2] -> [5];
    // [2][0] -> [5];

    VecD<6> activeStress;
    activeStress[0] = factor * fiberDirection[0] * fiberDirection[0];
    activeStress[1] = factor * fiberDirection[1] * fiberDirection[1];
    activeStress[2] = factor * fiberDirection[2] * fiberDirection[2];
    activeStress[3] = factor * fiberDirection[0] * fiberDirection[1];
    activeStress[4] = factor * fiberDirection[1] * fiberDirection[2];
    activeStress[5] = factor * fiberDirection[0] * fiberDirection[2];

    LOG(DEBUG) << "dof " << dofNoLocal << ", lambda: " << lambda << ", lambdaRelative: " << lambdaRelative
      << ", pmax_: " << pmax_ << ", f: " << f << ", gamma: " << gamma << ", => factor: " << factor << ", fiberDirection: " << fiberDirection;

    // if lambda is not yet computed (before first computation), set active stress to zero
    if (fabs(lambda)  < 1e-12)
    {
      activeStress = VecD<6>{0.0};
    }

    LOG(DEBUG) << "set active stress to " << activeStress;
    activePK2StressVariable->setValue(dofNoLocal, activeStress, INSERT_VALUES);
  }

  activePK2StressVariable->zeroGhostBuffer();
  activePK2StressVariable->finishGhostManipulation();
  activePK2StressVariable->startGhostManipulation();
}

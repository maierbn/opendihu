#include "control/precice/surface_coupling/00_nested_solver.h"

#include <sstream>

namespace Control
{

template<typename T1, typename T2, typename T3>
std::shared_ptr<typename PreciceAdapterNestedSolver<Control::Coupling<Control::MultipleInstances<T1>,MuscleContractionSolver<T2,T3>>>::FunctionSpace>
PreciceAdapterNestedSolver<Control::Coupling<Control::MultipleInstances<T1>,MuscleContractionSolver<T2,T3>>>::
functionSpace(NestedSolverType &nestedSolver)
{
  return nestedSolver.timeStepping2().data().functionSpace();
}

template<typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<Control::MultipleInstances<T1>,MuscleContractionSolver<T2,T3>>>::
addDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                               std::vector<typename SpatialDiscretization::DirichletBoundaryConditionsBase<FunctionSpace,6>::ElementWithNodes> &dirichletBoundaryConditionElements)
{
  // add the dirichlet bc values
  bool overwriteBcOnSameDof = true;
  nestedSolver.timeStepping2().dynamicHyperelasticitySolver()->addDirichletBoundaryConditions(dirichletBoundaryConditionElements, overwriteBcOnSameDof);
}

//! update existing boundary conditions with new values
template<typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<Control::MultipleInstances<T1>,MuscleContractionSolver<T2,T3>>>::
updateDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                                  std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBoundaryConditionValues)
{
  nestedSolver.timeStepping2().dynamicHyperelasticitySolver()->updateDirichletBoundaryConditions(newDirichletBoundaryConditionValues);
}

template<typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<Control::MultipleInstances<T1>,MuscleContractionSolver<T2,T3>>>::
updateNeumannBoundaryConditions(NestedSolverType &nestedSolver,
                                std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,Quadrature::Gauss<3>,3>> neumannBoundaryConditions)
{
  nestedSolver.timeStepping2().dynamicHyperelasticitySolver()->hyperelasticitySolver().updateNeumannBoundaryConditions(neumannBoundaryConditions);
}

//! get the displacement and velocity vectors of the given local dof nos
template<typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<Control::MultipleInstances<T1>,MuscleContractionSolver<T2,T3>>>::
getDisplacementVelocityValues(NestedSolverType &nestedSolver, const std::vector<dof_no_t> &dofNosLocal,
                               std::vector<double> &displacementValues, std::vector<double> &velocityValues)
{
  // get the displacement values
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.timeStepping2().data().displacements()->getValues(dofNosLocal, values);

  // store displacement values in interleaved order (array of struct)
  int nVectors = values.size();
  displacementValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++)
  {
    displacementValues[3*i + 0] = values[i][0];
    displacementValues[3*i + 1] = values[i][1];
    displacementValues[3*i + 2] = values[i][2];
  }

  // get the velocity values
  values.clear();
  nestedSolver.timeStepping2().data().velocities()->getValues(dofNosLocal, values);

  // store velocity values in interleaved order (array of struct)
  nVectors = values.size();
  velocityValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++)
  {
    velocityValues[3*i + 0] = values[i][0];
    velocityValues[3*i + 1] = values[i][1];
    velocityValues[3*i + 2] = values[i][2];
  }
}

//! get the traction vectors of the given local dof nos
template<typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<Control::MultipleInstances<T1>,MuscleContractionSolver<T2,T3>>>::
getTractionValues(NestedSolverType &nestedSolver, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &tractionValues)
{
  /*std::vector<Vec3> values0;
  nestedSolver.timeStepping2().data().materialTraction()->getValuesWithoutGhosts(values0);

  std::stringstream s;
  for (int i = 0; i < values0.size(); i++)
  {
    s << " " << values0[i][2];
  }
  LOG(INFO) << values0.size() << " local traction values in total (MuscleContractionSolver), "
    << "\ndofNosLocal: " << dofNosLocal << "\nall values: " << values0 << "\nz values: " << s.str();
*/
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.timeStepping2().data().materialTraction()->getValues(dofNosLocal, values);

  int nVectors = values.size();
  tractionValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++)
  {
    tractionValues[3*i + 0] = values[i][0];
    tractionValues[3*i + 1] = values[i][1];
    tractionValues[3*i + 2] = values[i][2];
  }
}

template<typename T1, typename T2, typename T3>
Vec PreciceAdapterNestedSolver<Control::Coupling<Control::MultipleInstances<T1>,MuscleContractionSolver<T2,T3>>>::
currentState(NestedSolverType &nestedSolver)
{
  return nestedSolver.timeStepping2().dynamicHyperelasticitySolver()->currentState();
}

// --------------------------------------------------

template<typename Material>
std::shared_ptr<typename TimeSteppingScheme::DynamicHyperelasticitySolver<Material>::FunctionSpace>
PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
functionSpace(NestedSolverType &nestedSolver)
{
  return nestedSolver.data().functionSpace();
}

template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
addDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                               std::vector<typename SpatialDiscretization::DirichletBoundaryConditionsBase<FunctionSpace,6>::ElementWithNodes> &dirichletBoundaryConditionElements)
{
  // add the dirichlet bc values
  bool overwriteBcOnSameDof = true;
  nestedSolver.addDirichletBoundaryConditions(dirichletBoundaryConditionElements, overwriteBcOnSameDof);
}

//! update existing boundary conditions with new values
template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
updateDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                                  std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBoundaryConditionValues)
{
  nestedSolver.updateDirichletBoundaryConditions(newDirichletBoundaryConditionValues);
}

template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
updateNeumannBoundaryConditions(NestedSolverType &nestedSolver,
                                std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,Quadrature::Gauss<3>,3>> neumannBoundaryConditions)
{
  nestedSolver.hyperelasticitySolver().updateNeumannBoundaryConditions(neumannBoundaryConditions);
}

//! get the displacement and velocity vectors of the given local dof nos
template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
getDisplacementVelocityValues(NestedSolverType &nestedSolver, const std::vector<dof_no_t> &dofNosLocal,
                              std::vector<double> &displacementValues, std::vector<double> &velocityValues)
{
  // get the displacement values
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.data().displacements()->getValues(dofNosLocal, values);

  // store displacement values in interleaved order (array of struct)
  int nVectors = values.size();
  displacementValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++)
  {
    displacementValues[3*i + 0] = values[i][0];
    displacementValues[3*i + 1] = values[i][1];
    displacementValues[3*i + 2] = values[i][2];
  }

  // get the velocity values
  values.clear();
  nestedSolver.data().velocities()->getValues(dofNosLocal, values);

  // store velocity values in interleaved order (array of struct)
  nVectors = values.size();
  velocityValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++)
  {
    velocityValues[3*i + 0] = values[i][0];
    velocityValues[3*i + 1] = values[i][1];
    velocityValues[3*i + 2] = values[i][2];
  }
}

//! get the traction vectors of the given local dof nos
template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
getTractionValues(NestedSolverType &nestedSolver, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &tractionValues)
{
  /*std::vector<Vec3> values0;
  nestedSolver.hyperelasticitySolver().data().materialTraction()->getValuesWithoutGhosts(values0);

  std::stringstream s;
  for (int i = 0; i < values0.size(); i++)
  {
    s << " " << values0[i][2];
  }
  LOG(INFO) << values0.size() << " local traction values in total (DynamicHyperelasticitySolver), "
    << "\ndofNosLocal: " << dofNosLocal << "\nall values: " << values0 << "\nz values: " << s.str();
*/
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.hyperelasticitySolver().data().materialTraction()->getValues(dofNosLocal, values);

  int nVectors = values.size();
  tractionValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++)
  {
    tractionValues[3*i + 0] = values[i][0];
    tractionValues[3*i + 1] = values[i][1];
    tractionValues[3*i + 2] = values[i][2];
  }
}

template<typename Material>
Vec PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
currentState(NestedSolverType &nestedSolver)
{
  return nestedSolver.currentState();
}


}  // namespace

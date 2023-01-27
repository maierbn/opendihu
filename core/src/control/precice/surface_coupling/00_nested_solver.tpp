#include "control/precice/surface_coupling/00_nested_solver.h"

#include <sstream>

namespace Control
{

// --------------------------------------------------
// Coupling<T1,MuscleContractionSolver<T2,T3>>
template<typename T1, typename T2, typename T3>
std::shared_ptr<typename PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::FunctionSpace>
PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::
functionSpace(NestedSolverType &nestedSolver)
{
  return nestedSolver.timeStepping2().data().functionSpace();
}

template<typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::
addDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                               std::vector<typename SpatialDiscretization::DirichletBoundaryConditionsBase<FunctionSpace,6>::ElementWithNodes> &dirichletBoundaryConditionElements)
{
  // add the dirichlet bc values
  bool overwriteBcOnSameDof = true;
  nestedSolver.timeStepping2().dynamicHyperelasticitySolver()->addDirichletBoundaryConditions(dirichletBoundaryConditionElements, overwriteBcOnSameDof);
}

//! update existing boundary conditions with new values
template<typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::
updateDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                                  std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBoundaryConditionValues)
{
  nestedSolver.timeStepping2().dynamicHyperelasticitySolver()->updateDirichletBoundaryConditions(newDirichletBoundaryConditionValues);
}

template<typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::
updateNeumannBoundaryConditions(NestedSolverType &nestedSolver,
                                std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,Quadrature::Gauss<3>,3>> neumannBoundaryConditions)
{
  nestedSolver.timeStepping2().dynamicHyperelasticitySolver()->hyperelasticitySolver().updateNeumannBoundaryConditions(neumannBoundaryConditions);
}

//! get the displacement and velocity vectors of the given local dof nos
template<typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::
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
void PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::
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
Vec PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::
currentState(NestedSolverType &nestedSolver)
{
  //Vec dynamicHyperelasticSolverState;
  // std::vector<double> fastMonodomainSolverState;
  // fastMonodomainSolverState = nestedSolver.timeStepping1().currentState();
  return nestedSolver.timeStepping2().dynamicHyperelasticitySolver()->currentState();
}

template<typename T1, typename T2, typename T3>
std::shared_ptr<FieldVariable::FieldVariable<typename PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::FunctionSpace,9>>
PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::
deformationGradientField(NestedSolverType &nestedSolver)
{
  return nestedSolver.timeStepping2().dynamicHyperelasticitySolver()->hyperelasticitySolver().data().deformationGradient();
}

template<typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<T1,MuscleContractionSolver<T2,T3>>>::
reset(NestedSolverType &nestedSolver)
{
  //nestedSolver.timeStepping1().reset();
  nestedSolver.timeStepping2().reset();
}

// --------------------------------------------------
// DynamicHyperelasticitySolver

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
  LOG(INFO) << "add dirichlet BC \n";

  bool overwriteBcOnSameDof = true;
  nestedSolver.addDirichletBoundaryConditions(dirichletBoundaryConditionElements, overwriteBcOnSameDof);
}

//! update existing boundary conditions with new values
template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
updateDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                                  std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBoundaryConditionValues)
{
    LOG(INFO) << "update dirichlet BC \n";
  nestedSolver.updateDirichletBoundaryConditions(newDirichletBoundaryConditionValues);
}

template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
updateNeumannBoundaryConditions(NestedSolverType &nestedSolver,
                                std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,Quadrature::Gauss<3>,3>> neumannBoundaryConditions)
{
    LOG(INFO) << "update neumann BC \n";
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

  // // get the velocity values
  // values.clear();
  // nestedSolver.data().velocities()->getValues(dofNosLocal, values);

  // // store velocity values in interleaved order (array of struct)
  // nVectors = values.size();
  // velocityValues.resize(nVectors * 3);

  // for (int i = 0; i < nVectors; i++)
  // {
  //   velocityValues[3*i + 0] = values[i][0];
  //   velocityValues[3*i + 1] = values[i][1];
  //   velocityValues[3*i + 2] = values[i][2];
  // }
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

template<typename Material>
std::shared_ptr<FieldVariable::FieldVariable<typename PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::FunctionSpace, 9>>
PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
deformationGradientField(NestedSolverType &nestedSolver)
{
  return nestedSolver.hyperelasticitySolver().data().deformationGradient();
}

template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
reset(NestedSolverType &nestedSolver)
{
  //nestedSolver.reset();
}

// --------------------------------------------------
// QHyperelasticitySolver

template<typename Material>
std::shared_ptr<typename TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>::FunctionSpace>
PreciceAdapterNestedSolver<TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>>::
functionSpace(NestedSolverType &nestedSolver)
{
  return nestedSolver.data().functionSpace();
}

template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>>::
addDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                               std::vector<typename SpatialDiscretization::DirichletBoundaryConditionsBase<FunctionSpace,6>::ElementWithNodes> &dirichletBoundaryConditionElements)
{
  // add the dirichlet bc values
  bool overwriteBcOnSameDof = true;
  nestedSolver.addDirichletBoundaryConditions(dirichletBoundaryConditionElements, overwriteBcOnSameDof);
}

//! update existing boundary conditions with new values
template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>>::
updateDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                                  std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBoundaryConditionValues)
{

  LOG(INFO) << "update dirichlet BC \n";
  
  nestedSolver.updateDirichletBoundaryConditions(newDirichletBoundaryConditionValues);
}

template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>>::
updateNeumannBoundaryConditions(NestedSolverType &nestedSolver,
                                std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,Quadrature::Gauss<3>,3>> neumannBoundaryConditions)
{
  LOG(INFO) << "update neumann BC \n";
  nestedSolver.hyperelasticitySolver().updateNeumannBoundaryConditions(neumannBoundaryConditions);
}

//! get the displacement and velocity vectors of the given local dof nos
template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>>::
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

  velocityValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++)
  {
    velocityValues[3*i + 0] = 0.0;
    velocityValues[3*i + 1] = 0.0;
    velocityValues[3*i + 2] = 0.0;
  }
}

//! get the traction vectors of the given local dof nos
template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>>::
getTractionValues(NestedSolverType &nestedSolver, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &tractionValues)
{
  /*std::vector<Vec3> values0;
  nestedSolver.hyperelasticitySolver().data().materialTraction()->getValuesWithoutGhosts(values0);

  std::stringstream s;
  for (int i = 0; i < values0.size(); i++)
  {
    s << " " << values0[i][2];
  }
  LOG(INFO) << values0.size() << " local traction values in total (QuasistaticHyperelasticitySolver), "
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
Vec PreciceAdapterNestedSolver<TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>>::
currentState(NestedSolverType &nestedSolver)
{
  return nestedSolver.currentState();
}

template<typename Material>
std::shared_ptr<FieldVariable::FieldVariable<typename PreciceAdapterNestedSolver<TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>>::FunctionSpace, 9>>
PreciceAdapterNestedSolver<TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>>::
deformationGradientField(NestedSolverType &nestedSolver)
{
  return nestedSolver.hyperelasticitySolver().data().deformationGradient();
}

template<typename Material>
void PreciceAdapterNestedSolver<TimeSteppingScheme::QuasistaticHyperelasticitySolver<Material>>::
reset(NestedSolverType &nestedSolver)
{
  nestedSolver.reset();
}


// --------------------------------------------------
// HyperelasticitySolver

template<typename Material>
std::shared_ptr<typename SpatialDiscretization::HyperelasticitySolver<Material>::FunctionSpace>
PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<Material>>::
functionSpace(NestedSolverType &nestedSolver)
{
  return nestedSolver.data().functionSpace();
}

template<typename Material>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<Material>>::
addDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                               std::vector<typename SpatialDiscretization::DirichletBoundaryConditionsBase<FunctionSpace,6>::ElementWithNodes> &dirichletBoundaryConditionElements)
{
  // transform the ElementWithNodes variable from one with 6 components (displacements + velocities) to one with only 3 components
  std::vector<typename SpatialDiscretization::DirichletBoundaryConditionsBase<FunctionSpace,3>::ElementWithNodes> dirichletBoundaryConditionElements3;
  for (typename SpatialDiscretization::DirichletBoundaryConditionsBase<FunctionSpace,6>::ElementWithNodes &element6 : dirichletBoundaryConditionElements)
  {
    typename SpatialDiscretization::DirichletBoundaryConditionsBase<FunctionSpace,3>::ElementWithNodes element3;

    // copy elementNo from element6 to element3
    element3.elementNoLocal = element6.elementNoLocal;

    // copy all bc entries from element6 to element3, use the first 3 components for every prescribed value
    for (std::map<int,std::array<double,6>>::const_iterator iter = element6.elementalDofIndex.cbegin(); iter != element6.elementalDofIndex.cend(); iter++)
    {
      int elementalDofIndex = iter->first;
      std::array<double,3> value = {
        iter->second[0],
        iter->second[1],
        iter->second[2]
      };
      element3.elementalDofIndex.insert(std::pair<int,std::array<double,3>>(elementalDofIndex,value));
    }

    dirichletBoundaryConditionElements3.push_back(element3);
  }

  // add the dirichlet bc values
  bool overwriteBcOnSameDof = true;
  nestedSolver.addDirichletBoundaryConditions(dirichletBoundaryConditionElements3, overwriteBcOnSameDof);
}

//! update existing boundary conditions with new values
template<typename Material>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<Material>>::
updateDirichletBoundaryConditions(NestedSolverType &nestedSolver,
                                  std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBoundaryConditionValues)
{
  // transform the ElementWithNodes variable from one with 6 components (displacements + velocities) to one with only 3 components
  std::vector<std::pair<global_no_t,std::array<double,3>>> newDirichletBoundaryConditionValues3;
  for (std::pair<global_no_t,std::array<double,6>> &element6 : newDirichletBoundaryConditionValues)
  {
    std::pair<global_no_t,std::array<double,3>> element3;

    element3.first = element6.first;

    element3.second[0] = element6.second[0];
    element3.second[1] = element6.second[1];
    element3.second[2] = element6.second[2];

    newDirichletBoundaryConditionValues3.push_back(element3);
  }

  nestedSolver.updateDirichletBoundaryConditions(newDirichletBoundaryConditionValues3);
}

template<typename Material>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<Material>>::
updateNeumannBoundaryConditions(NestedSolverType &nestedSolver,
                                std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,Quadrature::Gauss<3>,3>> neumannBoundaryConditions)
{
  nestedSolver.updateNeumannBoundaryConditions(neumannBoundaryConditions);
}

//! get the displacement and velocity vectors of the given local dof nos
template<typename Material>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<Material>>::
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
  velocityValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++)
  {
    velocityValues[3*i + 0] = 0.0;
    velocityValues[3*i + 1] = 0.0;
    velocityValues[3*i + 2] = 0.0;
  }
}

//! get the traction vectors of the given local dof nos
template<typename Material>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<Material>>::
getTractionValues(NestedSolverType &nestedSolver, const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &tractionValues)
{
  /*std::vector<Vec3> values0;
  nestedSolver.data().materialTraction()->getValuesWithoutGhosts(values0);

  std::stringstream s;
  for (int i = 0; i < values0.size(); i++)
  {
    s << " " << values0[i][2];
  }
  LOG(INFO) << values0.size() << " local traction values in total (), "
    << "\ndofNosLocal: " << dofNosLocal << "\nall values: " << values0 << "\nz values: " << s.str();
*/
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.data().materialTraction()->getValues(dofNosLocal, values);

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
Vec PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<Material>>::
currentState(NestedSolverType &nestedSolver)
{
  return nestedSolver.currentState();
}

template<typename Material>
std::shared_ptr<FieldVariable::FieldVariable<typename PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<Material>>::FunctionSpace, 9>>
PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<Material>>::
deformationGradientField(NestedSolverType &nestedSolver)
{
  return nestedSolver.data().deformationGradient();
}

template<typename Material>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<Material>>::
reset(NestedSolverType &nestedSolver)
{
  nestedSolver.reset();
}

}  // namespace

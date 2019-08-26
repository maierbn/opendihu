#include "spatial_discretization/mechanics_solver.h"

namespace SpatialDiscretization
{
template<typename FunctionSpaceTypeType>
MechanicsSolver<FunctionSpaceType>::
MechanicsSolver(DihuContext context) :
  context_(context["MechanicsSolver"]), data_(context_), neumannBoundaryConditions_(context_)
{
  // create function space from settings
  std::shared_ptr<FunctionSpaceType> functionSpace = context_.meshManager()->functionSpace<FunctionSpaceTypeType>(context.getPythonConfig());

  // store mesh in data
  data_.setFunctionSpaceType(functionSpace);
}

template<typename FunctionSpaceType>
void MechanicsSolver<FunctionSpaceType>::
initialize()
{
  data_.initialize();

  setStiffnessMatrix();
  setRightHandSide();

  initialized_ = true;
}

template<typename FunctionSpaceType>
void MechanicsSolver<FunctionSpaceType>::
run()
{

}

template<typename FunctionSpaceType>
void MechanicsSolver<FunctionSpaceType>::
setStiffnessMatrix()
{
}

template<typename FunctionSpaceType>
void MechanicsSolver<FunctionSpaceType>::
setRightHandSide()
{
  std::string boundaryConditionsConfigKey = "tractionBoundaryConditions";
  neumannBoundaryConditions_.initialize(context_.getPythonConfig(), data_.functionSpace(), boundaryConditionsConfigKey);

  // rhs = -∫ (T_L*phi_L)_a δu_aM phi_M ds,   a = 1,2,3, L,M ∈ surface dofs
  data_.setExternalForcesRightHandSide(neumannBoundaryConditions_.rhs());
}

template<typename FunctionSpaceType>
void MechanicsSolver<FunctionSpaceType>::reset()
{
  initialized_ = false;
}

template<typename FunctionSpaceType>
bool MechanicsSolver<FunctionSpaceType>::
knowsMeshType()
{
  return true;
}

template<typename FunctionSpaceType>
MechanicsSolver<FunctionSpaceType>::Data &MechanicsSolver<FunctionSpaceType>::
data()
{
  return data_;
}

template<typename FunctionSpaceType>
MechanicsSolver<FunctionSpaceType>::TransferableSolutionDataType MechanicsSolver<FunctionSpaceType>::
getSolutionForTransfer()
{
  return data_.getSolutionForTransfer();
}

}  // namespace


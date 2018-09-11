#include "time_stepping_scheme/time_stepping_scheme_ode.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"
#include "utility/boundary_conditions.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOde<DiscretizableInTimeType>::
initializeBoundaryConditions()
{
  if (PythonUtility::hasKey(this->specificSettings_, "dirichletBoundaryConditions"))
  {
    BoundaryConditions::parseBoundaryConditions<typename DiscretizableInTimeType::FunctionSpace>(
      this->specificSettings_, this->data_->functionSpace(), boundaryConditions_);

    // determine if the BC indices in the config are given for global or local dof nos
    bool inputMeshIsGlobal = PythonUtility::getOptionBool(this->specificSettings_, "inputMeshIsGlobal", true);

    const int D = DiscretizableInTimeType::FunctionSpace::dim();
    const int nDofsPerNode = DiscretizableInTimeType::FunctionSpace::nDofsPerNode();

    // split the vector of pair values in one vector for indices and one vector for values
    boundaryConditionDofs_.reserve(boundaryConditions_.size());
    boundaryConditionValues_.reserve(boundaryConditions_.size());
    for (int i = 0; i < boundaryConditions_.size(); i++)
    {
      dof_no_t boundaryConditionDof = boundaryConditions_[i].first;

      // transfer global dofs to local dofs
      if (inputMeshIsGlobal)
      {
        bool isLocalNonGhost = true;
        global_no_t nodeGlobalNo = boundaryConditions_[i].first / nDofsPerNode;
        std::array<int,D> localCoordinates = this->data_->functionSpace()->meshPartition()->getLocalCoordinates(nodeGlobalNo, isLocalNonGhost);

        if (!isLocalNonGhost)
          continue;

        boundaryConditionDof = this->data_->functionSpace()->getNodeNo(localCoordinates);
        VLOG(1) << "transform global dof no " << boundaryConditions_[i].first << " to local coordinates " << localCoordinates << ", then to local dof no " << boundaryConditionDof;
      }

      boundaryConditionDofs_.push_back(boundaryConditionDof);
      boundaryConditionValues_.push_back(boundaryConditions_[i].second);
    }

    VLOG(1) << "initialized boundary conditions for time stepping scheme: dofs: " << boundaryConditionDofs_ << ", values: " << boundaryConditionValues_;
  }
}

template<typename DiscretizableInTimeType>
void TimeSteppingSchemeOde<DiscretizableInTimeType>::
applyBoundaryConditions()
{
  //std::shared_ptr<Mesh::Mesh> mesh = discretizableInTime_.mesh();
  this->data_->solution().setValues(boundaryConditionDofs_, boundaryConditionValues_);

  VLOG(1) << "applied boundary conditions in time stepping scheme:  dofs: " << boundaryConditionDofs_ << ", values: " << boundaryConditionValues_;
}

} // namespace

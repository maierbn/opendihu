#include "data_management/specialized_solver/quasi_static_nonlinear_elasticity_chaste.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data
{

template<typename FunctionSpace>
QuasiStaticNonlinearElasticityChaste<FunctionSpace>::
QuasiStaticNonlinearElasticityChaste(DihuContext context) :
  Data<FunctionSpace>::Data(context)
{
}

template<typename FunctionSpace>
void QuasiStaticNonlinearElasticityChaste<FunctionSpace>::
initialize()
{
  // call initialize of base class
  Data<FunctionSpace>::initialize();
}

template<typename FunctionSpace>
void QuasiStaticNonlinearElasticityChaste<FunctionSpace>::
createPetscObjects()
{
  LOG(DEBUG) << "QuasiStaticNonlinearElasticityChaste::createPetscObjects";

  assert(this->functionSpace_);

  this->activation_ = this->functionSpace_->template createFieldVariable<1>("activation");
  this->activeStress_ = this->functionSpace_->template createFieldVariable<9>("activeStress");
  this->displacement_ = this->functionSpace_->template createFieldVariable<3>("displacement");
}

template<typename FunctionSpace>
std::shared_ptr<typename QuasiStaticNonlinearElasticityChaste<FunctionSpace>::FieldVariableType>
QuasiStaticNonlinearElasticityChaste<FunctionSpace>::
activation()
{
  return this->activation_;
}

template<typename FunctionSpace>
std::shared_ptr<typename QuasiStaticNonlinearElasticityChaste<FunctionSpace>::VectorFieldVariableType>
QuasiStaticNonlinearElasticityChaste<FunctionSpace>::
displacement()
{
  return this->displacement_;
}

template<typename FunctionSpace>
std::shared_ptr<typename QuasiStaticNonlinearElasticityChaste<FunctionSpace>::StressFieldVariableType>
QuasiStaticNonlinearElasticityChaste<FunctionSpace>::
activeStress()
{
  return this->activeStress_;
}

template<typename FunctionSpace>
void QuasiStaticNonlinearElasticityChaste<FunctionSpace>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << *this->activation_;
}

template<typename FunctionSpace>
typename QuasiStaticNonlinearElasticityChaste<FunctionSpace>::OutputFieldVariables QuasiStaticNonlinearElasticityChaste<FunctionSpace>::
getOutputFieldVariables()
{
  // these field variables will be written to output files
  return std::tuple_cat(
    std::tuple<std::shared_ptr<FieldVariableType>>(this->activation_),
    std::tuple<std::shared_ptr<StressFieldVariableType>>(this->activeStress_),
    std::tuple<std::shared_ptr<VectorFieldVariableType>>(this->displacement_)
  );
}

} // namespace

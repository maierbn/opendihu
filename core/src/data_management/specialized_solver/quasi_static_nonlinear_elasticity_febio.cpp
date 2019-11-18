#include "data_management/specialized_solver/quasi_static_nonlinear_elasticity_febio.h"

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data
{

QuasiStaticNonlinearElasticityFebio::
QuasiStaticNonlinearElasticityFebio(DihuContext context) :
  Data<FunctionSpace>::Data(context)
{
}

void QuasiStaticNonlinearElasticityFebio::
initialize()
{
  // call initialize of base class
  Data<FunctionSpace>::initialize();
}

void QuasiStaticNonlinearElasticityFebio::
createPetscObjects()
{
  LOG(DEBUG) << "QuasiStaticNonlinearElasticityFebio::createPetscObjects";

  assert(this->functionSpace_);

  this->activation_ = this->functionSpace_->template createFieldVariable<1>("activation");
}

std::shared_ptr<typename QuasiStaticNonlinearElasticityFebio::FieldVariableType>
QuasiStaticNonlinearElasticityFebio::
activation()
{
  return this->activation_;
}

void QuasiStaticNonlinearElasticityFebio::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << *this->activation_;
}

typename QuasiStaticNonlinearElasticityFebio::FieldVariablesForOutputWriter QuasiStaticNonlinearElasticityFebio::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  return std::tuple<std::shared_ptr<FieldVariableType>>(this->activation_);
}

} // namespace

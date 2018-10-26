#include "data_management/finite_element_method/finite_elements_mixed.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"
#include "function_space/function_space.h"
#include "mesh/unstructured_deformable.h"
#include "basis_function/hermite.h"

namespace Data
{

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename Term>
FieldVariable::FieldVariable<LowOrderFunctionSpaceType,1> &FiniteElements<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,Term>::
pressure()
{
  return *this->pressure_;
}


template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename Term>
typename FiniteElements<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,Term>::OutputFieldVariables FiniteElements<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,Term>::
getOutputFieldVariables()
{
  std::shared_ptr<FieldVariable::FieldVariable<HighOrderFunctionSpaceType,3>> actualGeometryField
    = std::make_shared<FieldVariable::FieldVariable<HighOrderFunctionSpaceType,3>>(this->functionSpace->geometryField());

  return OutputFieldVariables(
    this->geometryReference_,
    actualGeometryField,
    this->displacements_,
    this->pressure_,
    this->residual_,
    this->externalVirtualWork_
  );
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename Term>
void FiniteElements<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,Term>::
setFunctionSpace(std::shared_ptr<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>> mixedMesh)
{
  mixedMesh_ = mixedMesh;

  // store high order mesh as mesh_
  this->functionSpace = mixedMesh_->highOrderFunctionSpace();
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename Term>
std::shared_ptr<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>> FiniteElements<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,Term>::
mixedMesh()
{
  return mixedMesh_;
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename Term>
void FiniteElements<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,Term>::
initialize()
{
  // call initialize of parent class
  FiniteElementsSolidMechanics<HighOrderFunctionSpaceType,Term>::initialize();

  LOG(DEBUG) << "functionSpace has geometry field: " << this->functionSpace->hasGeometryField();
  initializeFieldVariables();
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename Term>
void FiniteElements<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,Term>::
initializeFieldVariables()
{
  std::vector<std::string> unnamedSingleComponent({"0"});
  this->pressure_ = std::make_shared<FieldVariable::FieldVariable<LowOrderFunctionSpaceType,1>>(this->mixedMesh_->lowOrderFunctionSpace(), "pressure", unnamedSingleComponent);
}

template<typename LowOrderFunctionSpaceType,typename HighOrderFunctionSpaceType,typename Term>
const dof_no_t FiniteElements<FunctionSpace::Mixed<LowOrderFunctionSpaceType,HighOrderFunctionSpaceType>,Term>::
getTangentStiffnessMatrixNRows()
{
  const int D = HighOrderFunctionSpaceType::dim();
  return this->mixedMesh_->highOrderFunctionSpace()->nDofsLocal() * D + this->mixedMesh_->lowOrderFunctionSpace()->nDofsLocal();
}

} // namespace Data

#include "data_management/specialized_solver/muscle_contraction_solver.h"

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

template<typename FunctionSpaceType>
MuscleContractionSolver<FunctionSpaceType>::
MuscleContractionSolver(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType>
void MuscleContractionSolver<FunctionSpaceType>::
initialize()
{
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();

  // create th output connector data object
  outputConnectorData_ = std::make_shared<OutputConnectorDataType>();

  // add all needed field variables to be transferred

  // add λ, λdot and γ as output connectors
  outputConnectorData_->addFieldVariable(this->lambda_);
  outputConnectorData_->addFieldVariable(this->lambdaDot_);
  outputConnectorData_->addFieldVariable(this->gamma_);

  // There is addFieldVariable(...) and addFieldVariable2(...) for the two different field variable types,
  // Refer to "data_management/output_connector_data.h" for details.
}

template<typename FunctionSpaceType>
void MuscleContractionSolver<FunctionSpaceType>::
createPetscObjects()
{
  assert(this->functionSpace_);

  // Here, the actual field variables will be created.
  // The string is the name of the field variable. It will also be used in the VTK output files.
  this->gamma_     = this->functionSpace_->template createFieldVariable<1>("γ");
  this->lambda_    = this->functionSpace_->template createFieldVariable<1>("λ");
  this->lambdaDot_ = this->functionSpace_->template createFieldVariable<1>("λdot");
}

template<typename FunctionSpaceType>
void MuscleContractionSolver<FunctionSpaceType>::
setFieldVariables(std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType> displacements,
                  std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType> velocities,
                  std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::StressFieldVariableType> activePK2Stress,
                  std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::StressFieldVariableType> pK2Stress,
                  std::shared_ptr<MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType> fiberDirection,
                  bool setGeometryFieldForTransfer)
{
  displacements_ = displacements;
  velocities_ = velocities;
  activePK2Stress_ = activePK2Stress;
  pK2Stress_ = pK2Stress;
  fiberDirection_ = fiberDirection;

  if (setGeometryFieldForTransfer)
  {
    outputConnectorData_->addGeometryField(std::make_shared<typename FunctionSpaceType::GeometryFieldType>(this->displacements_->functionSpace()->geometryField()));
  }
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> MuscleContractionSolver<FunctionSpaceType>::
lambda()
{
  return this->lambda_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> MuscleContractionSolver<FunctionSpaceType>::
lambdaDot()
{
  return this->lambdaDot_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> MuscleContractionSolver<FunctionSpaceType>::
gamma()
{
  return this->gamma_;
}

template<typename FunctionSpaceType>
std::shared_ptr<typename MuscleContractionSolver<FunctionSpaceType>::OutputConnectorDataType> MuscleContractionSolver<FunctionSpaceType>::
getOutputConnectorData()
{
  // return the output connector data object
  return this->outputConnectorData_;
}

template<typename FunctionSpaceType>
typename MuscleContractionSolver<FunctionSpaceType>::FieldVariablesForOutputWriter MuscleContractionSolver<FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files by the output writer

  // get the geometry field, which is always needed, from the function space
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  return std::make_tuple(
    geometryField,
    this->lambda_,           //< relative fiber stretch
    this->lambdaDot_,        //< contraction velocity
    this->gamma_,            //< gamma, the homogenized stress
    this->displacements_,    //< u, the displacements
    this->velocities_,       //< v, the velocities
    this->activePK2Stress_,  //< the symmetric PK2 stress tensor of the active contribution in Voigt notation
    this->pK2Stress_,        //< the symmetric PK2 stress tensor in Voigt notation
    this->fiberDirection_    //< direction of fibers at current point

  );
}

} // namespace

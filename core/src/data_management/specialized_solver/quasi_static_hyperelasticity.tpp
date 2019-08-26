#include "data_management/specialized_solver/quasi_static_hyperelasticity.h"

namespace Data
{

template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
  QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
QuasiStaticHyperelasticity(DihuContext context) :
  Data<DisplacementsFunctionSpace>::Data(context)
{
}

template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
initialize()
{
  // call initialize of base class, this calls createPetscObjects
  Data<DisplacementsFunctionSpace>::initialize();
}

template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
createPetscObjects()
{
  LOG(DEBUG) << "QuasiStaticHyperelasticity::createPetscObject";

  assert(this->displacementsFunctionSpace_);
  assert(this->pressureFunctionSpace_);
  assert(this->functionSpace_);


  std::vector<std::string> displacementsComponentNames({"x","y","z"});
  displacements_           = this->displacementsFunctionSpace_->template createFieldVariable<3>("u", displacementsComponentNames);     //< u, the displacements
  displacementsLinearMesh_ = this->pressureFunctionSpace_->template createFieldVariable<3>("uLin", displacementsComponentNames);     //< u, the displacements
  pressure_                = this->pressureFunctionSpace_->template createFieldVariable<1>("p");     //<  p, the pressure variable

  std::vector<std::string> componentNames({"S_11", "S_22", "S_33", "S_12", "S_13", "S_23"});
  pK2Stress_               = this->displacementsFunctionSpace_->template createFieldVariable<6>("PK2-Stress (Voigt)", componentNames);     //<  the symmetric PK2 stress tensor in Voigt notation
}


//! field variable of geometryReference_
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::DisplacementsFieldVariableType> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
geometryReference()
{
  return this->geometryReference_;
}

//! field variable of u
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::DisplacementsFieldVariableType> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
displacements()
{
  return this->displacements_;
}

//! field variable of p
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::PressureFieldVariableType> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
pressure()
{
  return this->pressure_;
}

//! field variable of S
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::StressFieldVariableType> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
pK2Stress()
{
  return this->pK2Stress_;
}

template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
updateGeometry()
{
  PetscErrorCode ierr;

  // update quadratic function space geometry
  // w = alpha * x + y, VecWAXPY(w, alpha, x, y)
  ierr = VecWAXPY(this->displacementsFunctionSpace_->geometryField().valuesGlobal(),
                  1, this->displacements_->valuesGlobal(), this->geometryReference_->valuesGlobal()); CHKERRV(ierr);

  // for displacements extract linear mesh from quadratic mesh
  std::vector<Vec3> displacementValues;
  this->displacements_->getValuesWithGhosts(displacementValues);

  node_no_t nNodesLocal[3] = {
    this->displacementsFunctionSpace_->meshPartition()->nNodesLocalWithGhosts(0),
    this->displacementsFunctionSpace_->meshPartition()->nNodesLocalWithGhosts(1),
    this->displacementsFunctionSpace_->meshPartition()->nNodesLocalWithGhosts(2)
  };

  std::vector<Vec3> linearMeshDisplacementValues(this->pressureFunctionSpace_->meshPartition()->nNodesLocalWithGhosts());
  int linearMeshIndex = 0;

  // loop over linear nodes in the quadratic mesh
  for (int k = 0; k < nNodesLocal[2]; k+=2)
  {
    for (int j = 0; j < nNodesLocal[1]; j+=2)
    {
      for (int i = 0; i < nNodesLocal[0]; i+=2, linearMeshIndex++)
      {
        int index = k*nNodesLocal[0]*nNodesLocal[1] + j*nNodesLocal[0] + i;

        linearMeshDisplacementValues[linearMeshIndex] = displacementValues[index];
      }
    }
  }

  displacementsLinearMesh_->setValuesWithGhosts(linearMeshDisplacementValues, INSERT_VALUES);

  // update linear function space geometry
  // w = alpha * x + y, VecWAXPY(w, alpha, x, y)
  ierr = VecWAXPY(this->pressureFunctionSpace_->geometryField().valuesGlobal(),
                  1, this->displacementsLinearMesh_->valuesGlobal(), this->geometryReferenceLinearMesh_->valuesGlobal()); CHKERRV(ierr);
}

//! set the function space object that discretizes the pressure field variable
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
setPressureFunctionSpace(std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace)
{
  pressureFunctionSpace_ = pressureFunctionSpace;

  // set the geometry field of the reference configuration as copy of the geometry field of the function space
  geometryReferenceLinearMesh_ = std::make_shared<DisplacementsLinearFieldVariableType>(pressureFunctionSpace_->geometryField(), "geometryReferenceLinearMesh");
  geometryReferenceLinearMesh_->setValues(pressureFunctionSpace_->geometryField());
}

//! set the function space object that discretizes the displacements field variable
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
setDisplacementsFunctionSpace(std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace)
{
  displacementsFunctionSpace_ = displacementsFunctionSpace;

  // also set the functionSpace_ variable which is from the parent class Data
  this->functionSpace_ = displacementsFunctionSpace;

  LOG(DEBUG) << "create geometry Reference";

  // set the geometry field of the reference configuration as copy of the geometry field of the function space
  geometryReference_ = std::make_shared<DisplacementsFieldVariableType>(displacementsFunctionSpace_->geometryField(), "geometryReference");
  geometryReference_->setValues(displacementsFunctionSpace_->geometryField());
}

//! get the displacements function space
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<DisplacementsFunctionSpace> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
displacementsFunctionSpace()
{
  return displacementsFunctionSpace_;
}

//! get the pressure function space
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<PressureFunctionSpace> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
pressureFunctionSpace()
{
  return pressureFunctionSpace_;
}

template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
print()
{
}

template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::OutputFieldVariables QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
getOutputFieldVariables()
{
  // these field variables will be written to output files
  return std::tuple_cat(
    std::tuple<std::shared_ptr<DisplacementsFieldVariableType>>(std::make_shared<typename DisplacementsFunctionSpace::GeometryFieldType>(this->displacementsFunctionSpace_->geometryField())), // geometry
    std::tuple<std::shared_ptr<DisplacementsFieldVariableType>>(this->displacements_),              // displacements_
    std::tuple<std::shared_ptr<StressFieldVariableType>>(this->pK2Stress_)         // pK2Stress_
  );
}

} // namespace

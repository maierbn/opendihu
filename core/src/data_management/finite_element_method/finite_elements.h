#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <tuple>

#include "data_management/data.h"
#include "data_management/finite_element_method/finite_elements_base.h"
#include "data_management/finite_element_method/diffusion_tensor_directional.h"
#include "data_management/finite_element_method/diffusion_tensor_constant.h"
#include "equation/diffusion.h"
#include "equation/linear_elasticity.h"
#include "equation/type_traits.h"
#include "data_management/finite_element_method/linear_stiffness.h"

namespace Data
{

/*
 * Data class for general finite elements method
 */
template<typename FunctionSpaceType,int nComponents,typename Term,typename = Term,typename = typename FunctionSpaceType::BasisFunction>
class FiniteElements :
  public FiniteElementsBase<FunctionSpaceType,nComponents>
{
public:

  //! constructor
  using FiniteElementsBase<FunctionSpaceType,nComponents>::FiniteElementsBase;

  //! initialize base class and diffusion tensor
  using FiniteElementsBase<FunctionSpaceType,nComponents>::initialize;
};

/*
 * partial specialization for term with constant diffusion tensor, e.g.
 * Equation::Dynamic::AnisotropicDiffusion or
 * Equation::Static::GeneralizedLaplace
 */
template<typename FunctionSpaceType,int nComponents,typename Term>
class FiniteElements<FunctionSpaceType,nComponents,Term,Equation::hasGeneralizedLaplaceOperator<Term>> :
  public FiniteElementsBase<FunctionSpaceType,nComponents>,
  public DiffusionTensorConstant<FunctionSpaceType>
{
public:

  //! constructor
  FiniteElements(DihuContext context);

  // !intialize base class and diffusion tensor
  virtual void initialize();
};


/** for directional diffusion use the diffusion tensor that depends upon a direction field
 */
template<typename FunctionSpaceType,int nComponents>
class FiniteElements<FunctionSpaceType,nComponents,Equation::Dynamic::DirectionalDiffusion> :
  public FiniteElementsBase<FunctionSpaceType,nComponents>,
  public DiffusionTensorDirectional<FunctionSpaceType>
{
public:
  //! constructor
  FiniteElements(DihuContext context);

  //! dummy method
  virtual void initialize(){};

  // !intialize base class and diffusion tensor which needs the direction field and the number of compartments in the multidomain context
  virtual void initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction,
                          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> spatiallyVaryingPrefactor,
                          bool useAdditionalDiffusionTensor = false);
};

/** for linear elasticity use the class that holds linear elastcity parameters, K and μ
 */
template<typename FunctionSpaceType,int nComponents,typename Term>
class FiniteElements<FunctionSpaceType,nComponents,Term,Equation::isLinearElasticity<Term>> :
  public LinearStiffness<FunctionSpaceType,nComponents>
{
public:
  //! constructor
  using LinearStiffness<FunctionSpaceType,nComponents>::LinearStiffness;

  typedef FunctionSpaceType FunctionSpace;
  using LinearStiffness<FunctionSpaceType,nComponents>::FieldVariablesForOutputWriter;

  //! initialize, store the reference geometry as copy of the current geometry
  void initialize();

  //! update the geometry of the mesh and function space with the displacements, scaled by scalingFactor
  void updateGeometry(double scalingFactor=1.0);

  //! compute the strain from the current displacement (which is the solution field variable)
  void computeStrain(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents*nComponents>> strain);

protected:

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> referenceGeometry_;  //< reference geometry
};

}  // namespace

#include "data_management/finite_element_method/finite_elements.tpp"

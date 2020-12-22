#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <memory>

#include "spatial_discretization/neumann_boundary_conditions/00_neumann_boundary_conditions_base.h"
#include "data_management/neumann_boundary_conditions.h"

namespace SpatialDiscretization
{

/**
 *  Neumann boundary conditions, that are specified for faces of elements. There are scalar boundary conditions representing normal flux.
 *  For solid mechanics, Neumann boundary conditions represent surface traction and can have a direction.
 *  Applying Neumann boundary conditions is simply substracting a rhs vector from the normal rhs.
 */

/** initializeRhs for mesh dimension 2 or 3
 */
template<typename FunctionSpaceType, typename QuadratureType, int nComponents>
class NeumannBoundaryConditionsInitializeRhs :
  public NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>
{
public:
  using NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::NeumannBoundaryConditionsBase;

  // initialize the rhs in data, this is called by initialize() or can be called manually to update the rhs
  virtual void initializeRhs();
};

/** Neumann boundary condition class for mesh dimension 2 or 3, nComponents > 1
 *  for traction BC, not flux
 *
 */
template<typename FunctionSpaceType, typename QuadratureType, int nComponents, typename DummyForTraits=typename FunctionSpaceType::Mesh>
class NeumannBoundaryConditions :
public NeumannBoundaryConditionsInitializeRhs<FunctionSpaceType,QuadratureType,nComponents>
{
public:
  using NeumannBoundaryConditionsInitializeRhs<FunctionSpaceType,QuadratureType,nComponents>::NeumannBoundaryConditionsInitializeRhs;

  typedef typename NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::ElementWithFaces ElementWithFaces;

  //! assign the deformation gradient field variable to the data object. It is needed for converting traction bc in current configuration to reference configuration
  void setDeformationGradientField(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,9>> deformationGradientField);

protected:

  //! parse an object of type ElementWithFaces from python config,
  //! example values:  {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
  //! if elementNoLocal is != -1, it will be used as value for the local element no, otherwise the value is parsed from config
  virtual ElementWithFaces parseElementWithFaces(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace, element_no_t elementNoLocal) override;

};

/** Neumann boundary condition class for mesh dimension 2 or 3, nComponents == 1
 *  for flux BC, not traction
 *
 */
template<typename FunctionSpaceType, typename QuadratureType>
class NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,1,Mesh::isNotDim<1,typename FunctionSpaceType::Mesh>> :
public NeumannBoundaryConditionsInitializeRhs<FunctionSpaceType,QuadratureType,1>
{
public:
  using NeumannBoundaryConditionsInitializeRhs<FunctionSpaceType,QuadratureType,1>::NeumannBoundaryConditionsInitializeRhs;

  typedef typename NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,1>::ElementWithFaces ElementWithFaces;

protected:

  //! parse an object of type ElementWithFaces from python config,
  //! example values:  {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
  //! if elementNoLocal is != -1, it will be used as value for the local element no, otherwise the value is parsed from config
  virtual ElementWithFaces parseElementWithFaces(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace, element_no_t elementNoLocal) override;

};

/** Neumann boundary condition class for 1D meshes
 * nComponents == FunctionSpaceType::dim() for traction boundary condition or nComponents = 1 for flux BC
 *
 */
template<typename FunctionSpaceType, typename QuadratureType, int nComponents>
class NeumannBoundaryConditions<FunctionSpaceType, QuadratureType, nComponents, Mesh::isDim<1,typename FunctionSpaceType::Mesh>> :
public NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>
{
public:
  using NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::NeumannBoundaryConditionsBase;

  typedef typename NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::ElementWithFaces ElementWithFaces;

  // initialize the rhs in data, this is called by initialize()
  virtual void initializeRhs();

protected:

  //! parse an object of type ElementWithFaces from python config,
  //! example values:  {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
  //! if elementNoLocal is != -1, it will be used as value for the local element no, otherwise the value is parsed from config
  virtual ElementWithFaces parseElementWithFaces(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace, element_no_t elementNoLocal) override;

};

template<typename FunctionSpaceType, typename QuadratureType, int nComponents, typename DummyForTraits>
std::ostream &operator<<(std::ostream &stream, const NeumannBoundaryConditions<FunctionSpaceType,QuadratureType,nComponents,DummyForTraits> &rhs);

} // namespace

#include "spatial_discretization/neumann_boundary_conditions/01_neumann_boundary_conditions.tpp"

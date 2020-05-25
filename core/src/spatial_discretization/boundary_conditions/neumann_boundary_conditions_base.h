#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <memory>

#include "spatial_discretization/boundary_conditions/boundary_conditions_base.h"
#include "data_management/neumann_boundary_conditions.h"

namespace SpatialDiscretization
{

/** A class that handles Neumann type boundary conditions.
 *  nComponents == FunctionSpaceType::dim() for traction boundary condition or nComponents = 1 for flux BC
  */
template<typename FunctionSpaceType, typename QuadratureType, int nComponents>
class NeumannBoundaryConditionsBase :
  public BoundaryConditionsBase<FunctionSpaceType, nComponents>
{
public:

  struct ElementWithFaces;

  //! constructor
  NeumannBoundaryConditionsBase(DihuContext context);

  //! parse config and extract boundary conditions specified under the given key, store in boundaryConditionElements_
  void initialize(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace,
                  std::string boundaryConditionsConfigKey) override;

  //! initialize directly
  void initialize(std::shared_ptr<FunctionSpaceType> functionSpace, const std::vector<ElementWithFaces> &boundaryConditionElements);

  //! return the negative rhs with contributions from Neumann boundary conditions, i.e. int_∂Ω T*δu_aL*phi_L ds
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rhs();

  /** An element with faces and corresponding prescribed Neumann boundary condition values
   *
   */
  struct ElementWithFaces
  {
    element_no_t elementNoLocal;                                     //< the local no of the element

    Mesh::face_t face;                                               //< face on which the Neumann BC is applied
    std::vector<std::pair<dof_no_t, VecD<nComponents>>> dofVectors;  //< <surface-local dof no, value>, nComponents == FunctionSpaceType::dim() for traction boundary condition or nComponents = 1 for flux BC
    std::vector<dof_no_t> surfaceDofs;                               //< dof nos of the volume element that correspond to the face / surface. These are different from the dofs in dofsVector which are numbered for the surface only, surfaceDofs are in the numbering of the volume element.
    // note, for flux BC, dofVectors[i].second is a VecD<1>
  };

protected:

  //! initialize the rhs in data, this is called by initialize()
  virtual void initializeRhs() = 0;

  //! parse an object of type ElementWithFaces frorm python config,
  //! example values:  {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
  //! if elementNoLocal is != -1, it will be used as value for the local element no, otherwise the value is parsed from config
  virtual ElementWithFaces parseElementWithFaces(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace, element_no_t elementNoLocal) = 0;

  std::shared_ptr<FunctionSpaceType> functionSpace_;          //< the function space of the computational mesh (not the edges/faces) in which the Neumann bc are set
  std::vector<ElementWithFaces> boundaryConditionElements_;   //< elements with prescribed Neumman boundary condition values

  bool divideNeumannBoundaryConditionValuesByTotalArea_;      //< if the value in dofVectors is to be divided by the total area of the surface of all elements that have neumann bc

  Data::NeumannBoundaryConditions<FunctionSpaceType, nComponents> data_;    //< data object that contains the rhs vector with the BC contribution
};

} // namespace

#include "spatial_discretization/boundary_conditions/neumann_boundary_conditions_base.tpp"

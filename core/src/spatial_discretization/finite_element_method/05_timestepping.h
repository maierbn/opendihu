#pragma once

#include "spatial_discretization/finite_element_method/04_rhs.h"

#include "mesh/mesh.h"
#include "time_stepping_scheme/discretizable_in_time.h"

namespace SpatialDiscretization
{

/** class used for timestepping as for diffusion equation
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class FiniteElementMethodTimeStepping :
  public FiniteElementMethodBaseRhs<BasisOnMeshType, QuadratureType, Term>,
  public DiscretizableInTime
{
public:
  FiniteElementMethodTimeStepping(DihuContext context);

  //! return the compile-time constant number of variable components of the solution field variable
  static constexpr int nComponents();

  //! proceed time stepping by computing output = stiffnessMatrix*input, output back in strong form
  void evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime);

  //! initialize for use with timestepping
  void initialize() override;

  //! return true because the object has a specified mesh type
  bool knowsMeshType();

  //! return the mesh that is stored in the data class
  std::shared_ptr<Mesh::Mesh> mesh();

  typedef BasisOnMeshType BasisOnMesh;   ///< the BasisOnMesh type needed for time stepping scheme

  friend class StiffnessMatrixTester;    ///< a class used for testing
protected:

  //! do nothing, needed for initialize of base class that is overridden anyway
  void setRightHandSide(){};

  //! Compute from the rhs in weak formulation the rhs vector in strong formulation
  void recoverRightHandSideStrongForm(Vec &result);

  //! check if the matrix and vector number of entries are correct such that stiffnessMatrix can be multiplied to rhs
  void checkDimensions(Mat &stiffnessMatrix, Vec &rhs);
};

};  // namespace

#include "spatial_discretization/finite_element_method/05_timestepping.tpp"

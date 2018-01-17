#pragma once

#include "spatial_discretization/finite_element_method/04_rhs.h"

#include "mesh/mesh.h"
#include "time_stepping_scheme/discretizable_in_time.h"

namespace SpatialDiscretization
{

/** class used for timestepping as for diffusion equation
 */
template<typename BasisOnMeshType, typename IntegratorType, typename Term>
class FiniteElementMethodTimeStepping :
  public FiniteElementMethodBaseRhs<BasisOnMeshType, IntegratorType, Term>,
  public DiscretizableInTime
{
public:
  FiniteElementMethodTimeStepping(const DihuContext &context);
 
  //! proceed time stepping by computing output = stiffnessMatrix*input, output back in strong form
  void evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime);
  
  //! initialize for use with timestepping
  void initialize() override;
  
  //! return true because the object has a specified mesh type
  bool knowsMeshType();
  
  typedef BasisOnMeshType BasisOnMesh;   ///< the BasisOnMesh type needed for time stepping scheme
  
  friend class StiffnessMatrixTester;    ///< a class used for testing 
protected:
 
  //! do nothing, needed for initialize of base class that is overridden anyway
  void setRightHandSide(){};
  
  //! Extract from the rhs in weak formulation the rhs vector in strong formulation
  void recoverRightHandSide(Vec &result);
  
  double relativeTolerance_;      ///< relative tolerance for solver 
};

};  // namespace

#include "spatial_discretization/finite_element_method/05_timestepping.tpp"
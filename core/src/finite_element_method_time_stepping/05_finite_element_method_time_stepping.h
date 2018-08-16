#pragma once

#include "spatial_discretization/finite_element_method/04_rhs.h"

#include "mesh/mesh.h"
#include "discretizable_in_time/discretizable_in_time.h"

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
  void evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime);
  
  //! timestepping rhs of equation Au^(t+1)=Rhs^(t), used for the case (M/dt-K)u^(t+1)=M/dtu^(t)
  void evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime);
  
  //! inverse of the lumped mass matrix
  void setInvLumMassMatrix();
  
  //! if the matrix is already initialized
  bool invLumMassMatrixSet();
  
  //! precomputes the system matrix A=I-M^(-1)K from the inverse of the mass matrix M^(-1) and stiffness matrix K
  void preComputeSystemMatrix(Mat &systemMatrix);
  
  //! precomputes the system matrix A=M/dt-K for variant 1 of the implicit Euler
  void preComputeSystemMatrix1();
  
  //! solves the linear system of equations resulting from the Implicit Euler method time discretization
  void solveLinearSystem(Vec &input, Vec &output);

  //! initialize for use with timestepping
  void initialize();
  
  //! initialize for use with timestepping
  void initialize(double timeStepWidth); 
  
  //! set initial values and return true
  //bool setInitialValues(Vec& initialValues) override;

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
  
private:
  bool invLumMassMatrixSet_=false;
  //std::shared_ptr<KSP> ksp_;  //to be once created not for each implicit time step
  
};

};  // namespace

#include "finite_element_method_time_stepping/05_finite_element_method_time_stepping.tpp"
#include "finite_element_method_time_stepping/05_finite_element_method_time_stepping_explicit.tpp"
#include "finite_element_method_time_stepping/05_finite_element_method_time_stepping_implicit.tpp"

#pragma once

#include <Python.h>  // has to be the first included header

#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "specialized_solver/solid_mechanics/hyperelasticity/02_hyperelasticity_solver.h"
#include "data_management/specialized_solver/quasistatic_hyperelasticity_solver.h"

namespace TimeSteppingScheme
{

/** This is a solver for a dynamic nonlinear finite elasticity problem.
 *  It uses the hyperelasticity solver for the static computations.
 *  The dynamic problem includes the static problem, only the right hand side changes to account for inertia effects.
 *  Furthermore, an integration scheme for 2nd order ODEs is needed. This is given by leap frog integration.
 *
 */
template<typename Term = Equation::SolidMechanics::MooneyRivlinIncompressible3D, bool withLargeOutput=true, typename MeshType = Mesh::StructuredDeformableOfDimension<3>>
class QuasistaticHyperelasticitySolver :
  public TimeSteppingScheme
{
public:

  typedef SpatialDiscretization::HyperelasticitySolver<Term,withLargeOutput,MeshType,6> HyperelasticitySolverType;    // the hyperelasticity solver that solves the nonlinear problem, 6 non-pressure components (u and v)
  typedef typename HyperelasticitySolverType::DisplacementsFunctionSpace DisplacementsFunctionSpace;
  typedef DisplacementsFunctionSpace FunctionSpace;
  typedef typename HyperelasticitySolverType::PressureFunctionSpace PressureFunctionSpace;

  typedef PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,Term,6> VecHyperelasticity;
  typedef PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace,Term,6> MatHyperelasticity;

  typedef Data::QuasistaticHyperelasticitySolver<DisplacementsFunctionSpace> Data;

  //! constructor
  QuasistaticHyperelasticitySolver(DihuContext context);

  //! advance simulation by the given time span
  void advanceTimeSpan(bool withOutputWritersEnabled = true);

  //! initialize everything for the simulation
  void initialize();

  //! run the whole simulation, repeatedly calls advanceTimeSpan
  void run();

  //! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
  void callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement = 1);

  //! return the data object, with the call to this method the output writers get the data to create their output files
  Data &data();

  //! get a reference to the underlying HyperelasticitySolver which has the material formulation and the nonlinear solver
  HyperelasticitySolverType &hyperelasticitySolver();

  //! set new dirichlet boundary condition values for existing dofs
  void updateDirichletBoundaryConditions(std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBCValues);

  //! add new dirichlet bc's, it is also possible to set new dofs that were not prescribed beforehand
  //! this calls addBoundaryConditions() of the dirichletBoundaryConditions_ object
  //! @param overwriteBcOnSameDof if existing bc dofs that are also in the ones to set newly should be overwritten, else they are not touched
  void addDirichletBoundaryConditions(std::vector<typename SpatialDiscretization::DirichletBoundaryConditions<DisplacementsFunctionSpace,6>::ElementWithNodes> &boundaryConditionElements, bool overwriteBcOnSameDof);

  //! get the Petsc Vec of the current state (uvp vector), this is needed to save and restore checkpoints from the PreciceAdapter
  Vec currentState();

private:

  //! set initial values for u and v from settings
  void setInitialValues();

  //! call the callback function to update dirichlet boundary condition values
  void callUpdateDirichletBoundaryConditionsFunction(double t);

  //! call the callback function to update Neumann boundary condition values
  void callUpdateNeumannBoundaryConditionsFunction(double t);

  //! recreate the rhs of the neumann boundary conditions, needed if traction was specified in current configuration
  void updateNeumannBoundaryConditions();

  //! compute the total bearing forces and moments at the bottom (z-) and top (z+) of the domain
  void computeBearingForcesAndMoments(double currentTime);

  HyperelasticitySolverType hyperelasticitySolver_;  //< hyperelasticity solver that solver the static problem
  Data data_;

  double density_;   //< density rho, used for inertia

  std::shared_ptr<VecHyperelasticity> uvp_;     //< combined vector of u,v and p values
  Vec internalVirtualWork_;                     //< internal virtual work, computed by the hyperelasticity solver
  Vec accelerationTerm_;                        //< contribution to virtual work from acceleration, computed by the hyperelasticity solver
  Vec externalVirtualWorkDead_;                 //< external virtual work, computed by the hyperelasticity solver, dead load i.e. constant over time

  bool inputMeshIsGlobal_;                      //< value of the setting "inputMeshIsGlobal", if the new dirichletBC values are given in global or local numbering

  PyObject *pythonUpdateDirichletBoundaryConditionsFunction_;     //< the callback function that updates dirichlet boundary conditions
  int updateDirichletBoundaryConditionsFunctionCallInterval_;     //< the interval with which the function will be called
  int updateDirichletBoundaryConditionsFunctionCallCount_ = 0;    //< the counter of number of calls to the updateDirichletBoundaryConditionsFunction

  PyObject *pythonUpdateNeumannBoundaryConditionsFunction_;       //< the callback function that updates neumann boundary conditions
  int updateNeumannBoundaryConditionsFunctionCallInterval_;       //< the interval with which the function will be called
  int updateNeumannBoundaryConditionsFunctionCallCount_ = 0;      //< the counter of number of calls to the updateNeumannBoundaryConditionsFunction

  PyObject *pythonTotalForceFunction_;                            //< the callback function that gets the total bearing forces and moments
  int pythonTotalForceFunctionCallInterval_;                      //< the interval with which the function will be called
  int pythonTotalForceFunctionCallCount_ = 0;                     //< the counter of number of calls to the pythonTotalForceFunction_

  std::vector<std::tuple<element_no_t,bool>> bottomTopElements_;  //< (elementNoLocal,isAtTopOfDomain) elements that are used to integrate total forces and moments
  std::string totalForceLogFilename_;                             //< filename of the log file that will contain the total bearing forces at top and bottom elements
  bool isTractionInCurrentConfiguration_;                         //< if traction is given in current configuration, then it has to be transformed to reference configuration in every timestep
  bool isReferenceGeometryInitialized_;                           //< if the copy of the reference geometry in the hyperelasticitySolver_ has already been set in the first timestep
};

}  // namespace

#include "specialized_solver/solid_mechanics/quasistatic_hyperelasticity/quasistatic_hyperelasticity_solver.tpp"

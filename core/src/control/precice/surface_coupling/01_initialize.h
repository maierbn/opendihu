#pragma once

#include <Python.h>  // has to be the first included header

#include "control/precice/surface_coupling/00_nested_solver.h"

#ifdef HAVE_PRECICE
#include "precice/SolverInterface.hpp"
#endif

namespace Control
{

/** Generic Precice adapter, can be configured to either prescribe Neumann or Dirichlet boundary conditions.
 */
template<typename NestedSolver>
class PreciceAdapterInitialize :
  public PreciceAdapterNestedSolver<NestedSolver>
{
public:

  //! constructor, gets the DihuContext object which contains all python settings
  PreciceAdapterInitialize(DihuContext context);

  //! initialize the object
  void initialize();

protected:

#ifdef HAVE_PRECICE

  using PreciceAdapterNestedSolver<NestedSolver>::FunctionSpace;

  /** a coupling mesh of precice, this does not have to coincide with an opendihu mesh
   */
  struct PreciceMesh
  {
    std::string preciceMeshName;              //< name of the precice mesh as used in the precice config XML file
    int preciceMeshId;                        //< mesh ID of precice of the surface mesh that contains all nodes
    std::vector<int> preciceVertexIds;        //< the vertex ids in precice of the geometry values
    std::vector<dof_no_t> dofNosLocal;        //< the local dof nos in the 3D mesh of the surface mesh nodes
    std::vector<double> geometryValuesSurface;  //< the geometry values, i.e., node positions

    int nNodesLocal;                          //< local number of nodes on the surface where the tendon is coupled
    enum {
      face2Minus,
      face2Plus
    } face;                                   //< the face of the 3D mesh where the 2D coupling mesh is, 2- = bottom, 2+ = top
  };

  /** a precice coupling participant to which the current solver is coupled to
   */
  struct PreciceData
  {
    std::string displacementsName;            //< precice name of the displacements variable, if any
    std::string velocitiesName;               //< precice name of the velocities variable, if any
    std::string tractionName;                 //< precice name of the traction variable, if any

    int preciceDataIdDisplacements;           //< precice id of the displacements variable
    int preciceDataIdVelocities;              //< precice id of the velocities variable
    int preciceDataIdTraction;                //< precice id of the traction variable

    bool average;

    enum {
      ioRead,
      ioWrite
    } ioType;                                 //< if this variable is to be written or read to other participants over precice

    enum {
      bcTypeDirichlet,
      bcTypeNeumann
    } boundaryConditionType;                  //< the type of the boundary condition to set, if ioType == ioRead

    std::shared_ptr<PreciceMesh> preciceMesh; //< the coupling mesh, this is derived from the option preciceMeshName
  };

  //! parse the options in "preciceMeshes" and store in variable preciceMeshes_
  void initializePreciceMeshes();

  //! initialize all meshes in precice from the variable preciceMeshes_
  void setMeshesInPrecice();

  //! parse the options in "preciceData" and initialize all variables in precice, store in variable preciceData_
  void initializePreciceData();

  //! initialize Dirichlet boundary conditions at all dofs that will get some prescribed values during coupling
  void initializeDirichletBoundaryConditions();
#endif

  DihuContext context_;                       //< object that contains the python config for the current context and the global singletons meshManager and solverManager
  PythonConfig specificSettings_;             //< python object containing the value of the python config dict with corresponding key

  NestedSolver nestedSolver_;                 //< the nested solver that is controlled by this class

#ifdef HAVE_PRECICE
  std::string preciceParticipantName_;        //< name of the participant as given in the precice config
  std::vector<std::shared_ptr<PreciceMesh>> preciceMeshes_;   //< all coupling meshes
  std::vector<PreciceData> preciceData_;      //< all precice variables "data"

  std::shared_ptr<precice::SolverInterface> preciceSolverInterface_;  //< the precice solver interface that makes all preCICE functionality accessible

  std::shared_ptr<typename PreciceAdapterNestedSolver<NestedSolver>::FunctionSpace> functionSpace_;   //< the function space of the nested solver

#endif

  bool ownRankIsInvolved_;                    //< if the own rank has part of a coupling surface and is involved in the coupling
  bool couplingEnabled_;                      //< if the coupling is enabled, if not it can be used for debugging, without precice
  double maximumPreciceTimestepSize_;         //< maximum timestep size that precice will allow for the current time step
  double timeStepWidth_;                      //< timestep width of the solver
  int timeStepOutputInterval_;                //< interval in which to output current time
  double scalingFactor_;                      //< a factor to scale the exchanged data prior to communication
  bool outputOnlyConvergedTimeSteps_;         //< option if the output should be written only for converged timesteps
  double endTimeIfCouplingDisabled_;          //< the end time that is used if the coupling is not enabled

  bool initialized_;                          //< if initialize() was already called
};

}  // namespace

#include "control/precice/surface_coupling/01_initialize.tpp"

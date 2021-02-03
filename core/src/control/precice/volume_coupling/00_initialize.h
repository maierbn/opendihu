#pragma once

#include <Python.h>  // has to be the first included header

#include "control/precice/volume_coupling/00_initialize.h"

#ifdef HAVE_PRECICE
#include "precice/SolverInterface.hpp"
#endif

namespace Control
{

/** Generic Precice adapter, can be configured to either prescribe Neumann or Dirichlet boundary conditions.
 */
template<typename NestedSolver>
class PreciceAdapterVolumeCouplingInitialize :
  public Runnable
{
public:

  //! constructor, gets the DihuContext object which contains all python settings
  PreciceAdapterVolumeCouplingInitialize(DihuContext context);

  //! initialize the object
  void initialize();

protected:

#ifdef HAVE_PRECICE

  typedef typename NestedSolver::FunctionSpace FunctionSpace;

  /** a coupling mesh of precice, this does not have to coincide with an opendihu mesh
   */
  struct PreciceMesh
  {
    std::string preciceMeshName;              //< name of the precice mesh as used in the precice config XML file
    int preciceMeshId;                        //< mesh ID of precice of the surface mesh that contains all nodes
    std::vector<int> preciceVertexIds;        //< the vertex ids in precice of the geometry values
    std::vector<dof_no_t> dofNosLocal;        //< the local dof nos in the 3D mesh of the surface mesh nodes
    std::vector<double> geometryValues;       //< the geometry values, i.e., node positions
    int nNodesLocal;                          //< local number of nodes
  };

  /** one field variable to be transferred, called `data:vector` or `data:scalar` in the precice config
   */
  struct PreciceData
  {
    std::string preciceDataName;              //< precice name of the variable, if any
    int preciceDataId;                        //< precice id of the variable
    std::string slotName;                     //< slot name as given in config, this is used to determine slotNo
    int slotNo;                               //< slot no that corresponds to this field variable

    std::string opendihuMeshName;             //< opendihu mesh name that is used for the geometry initialization
    bool isGeometryField;                     //< if the corresponding field variable is a geometry field

    enum {
      ioRead,
      ioWrite
    } ioType;                                 //< if this variable is to be written or read to other participants over precice

    std::shared_ptr<PreciceMesh> preciceMesh; //< the coupling mesh, this is derived from the option preciceMeshName
  };

  //! initialize all meshes in precice from the variable preciceMeshes_
  void setMeshesInPrecice();

  //! parse the options in "preciceData" and initialize all variables in precice, store in variable preciceData_
  void initializePreciceData();
#endif

  DihuContext context_;                       //< object that contains the python config for the current context and the global singletons meshManager and solverManager
  PythonConfig specificSettings_;             //< python object containing the value of the python config dict with corresponding key

  NestedSolver nestedSolver_;                 //< the nested solver that is controlled by this class

#ifdef HAVE_PRECICE
  std::string preciceParticipantName_;        //< name of the participant as given in the precice config
  std::vector<std::shared_ptr<PreciceMesh>> preciceMeshes_;   //< all coupling meshes
  std::vector<PreciceData> preciceData_;      //< all precice variables "data"

  std::shared_ptr<precice::SolverInterface> preciceSolverInterface_;  //< the precice solver interface that makes all preCICE functionality accessible

#endif

  bool ownRankIsInvolved_;                    //< if the own rank has part of a coupling surface and is involved in the coupling
  bool couplingEnabled_;                      //< if the coupling is enabled, if not it can be used for debugging, without precice
  double maximumPreciceTimestepSize_;         //< maximum timestep size that precice will allow for the current time step
  double timeStepWidth_;                      //< timestep width of the solver
  int timeStepOutputInterval_;                //< interval in which to output current time
  double scalingFactor_;                      //< a factor to scale the exchanged data prior to communication
  bool outputOnlyConvergedTimeSteps_;         //< option if the output should be written only for converged timesteps

  bool initialized_;                          //< if initialize() was already called
};

}  // namespace

#include "control/precice/volume_coupling/00_initialize.tpp"

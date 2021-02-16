#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "data_management/specialized_solver/muscle_contraction_solver.h"
#include "specialized_solver/solid_mechanics/dynamic_hyperelasticity/dynamic_hyperelasticity_solver.h"
#include "equation/mooney_rivlin_incompressible.h"
#include "data_management/specialized_solver/muscle_contraction_solver.h"

/** Solve the incompressible, transversely isotropic Mooney-Rivlin material with active stress contribution.
 *
 * This solver encapsulates the DynamicHyperelasticitySolver or HyperelasticitySolver (static). Which one to use can be chosen at runtime in the config.
 * This class adds functionality to compute and transfer active stresses, fiber stretches and contraction velocity and so on.
 */
template<typename MeshType=Mesh::StructuredDeformableOfDimension<3>,typename Term=Equation::SolidMechanics::TransverselyIsotropicMooneyRivlinIncompressibleActive3D,bool withLargeOutputFiles=true>
class MuscleContractionSolver :
  public Runnable,
  public TimeSteppingScheme::TimeSteppingScheme
{
public:
  typedef ::TimeSteppingScheme::DynamicHyperelasticitySolver<Term,withLargeOutputFiles,MeshType> DynamicHyperelasticitySolverType;
  typedef ::SpatialDiscretization::HyperelasticitySolver<Term,withLargeOutputFiles,MeshType> StaticHyperelasticitySolverType;

  //! make the DisplacementsFunctionSpace of the DynamicHyperelasticitySolver class available
  typedef typename DynamicHyperelasticitySolverType::DisplacementsFunctionSpace FunctionSpace;

  //! define the type of the data object,
  typedef typename Data::MuscleContractionSolver<FunctionSpace> Data;

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  //! Usually you define this type in the "Data" class and reuse it here.
  typedef typename Data::SlotConnectorDataType SlotConnectorDataType;

  //! constructor, gets the DihuContext object which contains all python settings
  MuscleContractionSolver(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_] (set by setTimeSpan(), take a look at time_stepping_scheme/00_time_stepping_scheme.h)
  void advanceTimeSpan(bool withOutputWritersEnabled = true);

  //! initialize time span from specificSettings_
  void initialize();

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary ("uninitialize")
  void reset();

  //! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
  void callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement = 1);

  //! return the data object of the timestepping scheme, with the call to this method the output writers get the data to create their output files
  Data &data();

  //! get a reference to the DynamicHyperelasticitySolverType
  std::shared_ptr<DynamicHyperelasticitySolverType> dynamicHyperelasticitySolver();

  //! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
  //! the transfer is done by the slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

protected:

  //! compute the active stress tensor field from lambda
  void computeActiveStress();

  //! compute λ and λ_dot for data transfer
  void computeLambda();

  //! update the geometry at the given meshes, map the own geometry field to the meshes
  void mapGeometryToGivenMeshes();

  //! create the mappings for the geometry field mapping between meshes
  void initializeMappingBetweenMeshes();

  std::shared_ptr<DynamicHyperelasticitySolverType> dynamicHyperelasticitySolver_;   //< the dynamic hyperelasticity solver that solves for the dynamic contraction
  std::shared_ptr<StaticHyperelasticitySolverType> staticHyperelasticitySolver_;     //< the static hyperelasticity solver that can be used for quasi-static solution

  Data data_;                                   //< the data object that holds all field variables
  OutputWriter::Manager outputWriterManager_;   //< manager object holding all output writers

  double pmax_;                                 //< settings of "Pmax" maximum active stress of the muscle
  bool isDynamic_;                              //< if the dynamic formulation or the quasi-static formulation is used
  bool enableForceLengthRelation_;              //< if the force-length relation factor f_l(λ_f) should be multiplied
  double lambdaDotScalingFactor_;               //< scaling factor for the computation of lambdaDot
  std::vector<std::string> meshNamesOfGeometryToMapTo_;   //< a list of mesh names which will get updated with the geometry

  bool initialized_;                            //< if initialize was already called
};

#include "specialized_solver/muscle_contraction_solver.tpp"
#include "specialized_solver/muscle_contraction_solver_compute.tpp"

#pragma once

#include <Python.h>  // has to be the first included header

namespace SpatialDiscretization
{

/** A solver for continuum mechanics problems with deformation
  */
template<typename FunctionSpaceType>
class MechanicsSolver :
  public Runnable
{
public:

  typedef Data::Mechanics<FunctionSpaceType> Data;

  //! constructor
  MechanicsSolver(DihuContext context);

  //! initialize components of the simulation
  void initialize();

  //! run the simulation
  void run();

  //! reset state
  void reset();

  //! return the data object
  Data &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  OutputConnectorDataType getOutputConnectorData();

protected:

  //! assemble the stiffness matrix in data for the solid mechanics problem
  void setStiffnessMatrix();

  //! assemble the right hand side of the system
  void setRightHandSide();

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  Data data_;  ///< the data object of the mechanics solver which stores all field variables and matrices

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writers

  NeumannBoundaryConditions<FunctionSpaceType,Quadrature::Gauss<3>,FunctionSpaceType::dim()> neumannBoundaryConditions_;   ///< neumann boundary conditions object, this parses the rhs with the external force

  bool initialized_;   ///< if this object was already initialized
  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
};

}  // namespace

#include "spatial_discretization/mechanics_solver.tpp"

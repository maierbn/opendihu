#pragma once

#include <Python.h>  // has to be the first included header

#include <vector>

#include "discretizable_in_time/discretizable_in_time.h"
#include "cellml/02_callback_handler.h"

/** The is a class that contains cellml equations and can be used with a time stepping scheme.
 *  The nStates template parameter specifies the number of state variables that should be used with the integrator.
 *  It is necessary that this value is fixed at compile time because the timestepping scheme needs to know which field variable types is has to construct.
 *  This class can also be computed easily in multiple instances along the nodes of a mesh. The number of instances is deduced from the mesh.
 * 
 *  The states values are not stored inside the class but in the time stepping scheme that is used to integrate the cellml problems.
 * 
 *  Naming:
 *   Intermediate (opendihu) = KNOWN (OpenCMISS) = Algebraic (OpenCOR)
 *   Parameter (opendihu, OpenCMISS) = KNOWN (OpenCMISS), in OpenCOR also algebraic
 *   Constant - these are constants that are only present in the source files
 *   State: state variable
 *   Rate: the time derivative of the state variable, i.e. the increment value in an explicit Euler stepping
 */
template <int nStates>
class CellmlAdapter :
  public CallbackHandler<nStates>
{
public:

  ///! constructor
  CellmlAdapter(DihuContext context);
  
  //! initialize callback functions and rhs
  void initialize();
  
  //! initialize timestepping
  void initialize(double timeStepWidth);
  
  //! reset the object to uninitialized state
  void reset();

  //! evaluate rhs
  void evaluateTimesteppingRightHandSideExplicit(Vec& input, Vec& output, int timeStepNo, double currentTime);
  
  //! evaluate rhs
  void evaluateTimesteppingRightHandSideImplicit(Vec& input, Vec& output, int timeStepNo, double currentTime);
  
  //! return false because the object is independent of mesh type
  bool knowsMeshType();

  //! return the mesh
  std::shared_ptr<Mesh::Mesh> mesh();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);
  
  //! set initial values and return true or don't do anything and return false
  bool setInitialValues(Vec &initialValues);
  
  //! get a vector with the names of the states
  void getComponentNames(std::vector<std::string> &stateNames) override;
  
 };

#include "cellml/03_cellml_adapter.tpp"

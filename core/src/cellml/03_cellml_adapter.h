#pragma once

#include <Python.h>  // has to be the first included header

#include <vector>

#include "time_stepping_scheme/discretizable_in_time.h"
#include "cellml/02_callback_handler.h"

/** The is a class that contains cellml equations and can be used with a time stepping scheme.
 *  The nStates template parameter specifies the number of state variables that should be used with the integrator.
 *  It is necessary that this value is fixed at compile time because the timestepping scheme needs to know which field variable types is has to construct.
 *  This class can also be computed easily in multiple instances along the nodes of a mesh.
 */
template <int nStates>
class CellmlAdapter :
  public CallbackHandler<nStates>,
  public DiscretizableInTime
{
public:

  ///! constructor
  CellmlAdapter(DihuContext context);
  
  //! initialize callback functions and rhs
  void initialize();
  
  //! reset state such that new initialization becomes necessary
  void reset();

  //! evaluate rhs
  void evaluateTimesteppingRightHandSide(Vec& input, Vec& output, int timeStepNo, double currentTime);
  
  //! return false because the object is independent of mesh type
  bool knowsMeshType();

  //! return the mesh
  std::shared_ptr<Mesh::Mesh> mesh();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);
  
  //! set initial values and return true or don't do anything and return false
  bool setInitialValues(Vec &initialValues);
  
private:

 };

#include "cellml/03_cellml_adapter.tpp"
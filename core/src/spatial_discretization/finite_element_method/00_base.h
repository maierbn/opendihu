#pragma once

#include "data_management/finite_elements.h"
#include "partition/rank_subset.h"

//#define QUADRATURE_TEST    ///< if evaluation of quadrature accuracy takes place
//#define EXACT_QUADRATURE Quadrature::Gauss<20>

namespace SpatialDiscretization
{

/**
 * Base class containing basic finite element functionality such as initializing and solving.
 * Further classes derive from this base class and add special functionality such as setting stiffness matrix, rhs and timestepping
 */
template<typename FunctionSpaceType,typename QuadratureType,typename Term>
class FiniteElementMethodBase : public SpatialDiscretization, public Runnable
{
public:
  FiniteElementMethodBase(DihuContext context);

  typedef ::Data::FiniteElements<FunctionSpaceType,Term> Data;
  typedef FunctionSpaceType FunctionSpace;

  // perform computation
  void run();

  //! initialize for use as laplace or poisson equation, not for timestepping
  virtual void initialize();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);
  
  //! get the stored mesh
  std::shared_ptr<Mesh::Mesh> mesh();

  //! get the data object
  Data &data();

  friend class StiffnessMatrixTester;    ///< a class used for testing
protected:

  //! read in rhs values from config and creates a FE rhs vector out of it
  virtual void setRightHandSide() = 0;

  //! setup stiffness matrix
  virtual void setStiffnessMatrix() = 0;
  
  //! setup mass matrix
  virtual void setMassMatrix() = 0;
  
  //! setup inverse of the lumped mass matrix
  //virtual void setInverseLumpedMassMatrix()=0;

  //! solve finite element linear system
  virtual void solve();

  //! after rhs is transferred to weak form this method is called and can be overriden later
  virtual void manipulateWeakRhs(){}

  //! modify the rhs to incorporate dirichlet boundary conditions
  virtual void applyBoundaryConditions() = 0;

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  Data data_;     ///< data object that holds all PETSc vectors and matrices
  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
};

};  // namespace

#include "spatial_discretization/finite_element_method/00_base.tpp"

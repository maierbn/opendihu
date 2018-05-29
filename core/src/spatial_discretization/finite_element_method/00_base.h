#pragma once

#include "data_management/finite_elements.h"

//#define QUADRATURE_TEST    ///< if evaluation of quadrature accuracy takes place
//#define EXACT_QUADRATURE Quadrature::Gauss<20>

namespace SpatialDiscretization
{

/**
 * Base class containing basic finite element functionality such as initializing and solving.
 * Further classes derive from this base class and add special functionality such as setting stiffness matrix, rhs and timestepping
 */
template<typename BasisOnMeshType,typename QuadratureType,typename Term>
class FiniteElementMethodBase : public SpatialDiscretization, public Runnable
{
public:
  FiniteElementMethodBase(DihuContext context);

  typedef ::Data::FiniteElements<BasisOnMeshType,Term> Data;
  typedef BasisOnMeshType BasisOnMesh;

  // perform computation
  void run();

  //! initialize for use as laplace or poisson equation, not for timestepping
  virtual void initialize();

  //! get the stored mesh
  std::shared_ptr<Mesh::Mesh> mesh();

  //! get the data object
  Data &data();

  friend class StiffnessMatrixTester;    ///< a class used for testing
protected:

  //! read in rhs values from config and creates a FE rhs vector out of it
  virtual void setRightHandSide() = 0;

  //! change the stiffness matrix such that Dirichlet boundary conditions are met, sets some rows/columns to 0 and the diagonal to 1, changes rhs accordingly
  virtual void applyBoundaryConditions();

  //! setup stiffness matrix
  virtual void setStiffnessMatrix() = 0;

  //! solve finite element linear system
  virtual void solve();

  //! after rhs is transferred to weak form this method is called and can be overriden later
  virtual void manipulateWeakRhs(){}

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  Data data_;     ///< data object that holds all PETSc vectors and matrices
  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
};

};  // namespace

#include "spatial_discretization/finite_element_method/00_base.tpp"

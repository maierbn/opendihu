#pragma once

#include "data_management/finite_elements.h"

namespace SpatialDiscretization
{
 

template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodBase : public SpatialDiscretization, public Runnable
{
public:
  FiniteElementMethodBase(const DihuContext &context);
  
  // perform computation
  void run();
  
  //! initialize for use as laplace or poisson equation, not for timestepping
  void initialize();
    
  //! get the stored mesh
  std::shared_ptr<Mesh::Mesh> mesh();
  
  friend class StiffnessMatrixTester;    ///< a class used for testing 
protected:
 
  virtual void setRightHandSide() = 0;
  virtual void applyBoundaryConditions();
  virtual void setStiffnessMatrix() = 0;
  virtual void solve();
  
  const DihuContext &context_;    ///< the context object containing everything to be stored
  Data::FiniteElements data_;     ///< data object that holds all PETSc vectors and matrices
  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
};
 

 
};  // namespace
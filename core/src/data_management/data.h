#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <vector>

#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "control/dihu_context.h"
#include "partition/rank_subset.h"

namespace Data
{

template<typename BasisOnMeshType>
class Data
{
public:

  typedef BasisOnMeshType BasisOnMesh;

  //! constructor
  Data(DihuContext context);

  //! destructor
  virtual ~Data();

  //! initialize the mesh with e.g. number of dimensions
  virtual void setMesh(std::shared_ptr<BasisOnMeshType> mesh);

  //! initialize, generate petsc objects, this has to be called after setMesh
  virtual void initialize();

  //! set the subset of ranks that will compute the work
  virtual void setRankSubset(Partition::RankSubset rankSubset);
  
  //! get the stored mesh
  const std::shared_ptr<BasisOnMeshType> mesh() const;

  //! return the total number of unknowns, i.e. degrees of freedom x number of components, this can be a multiple of the number of nodes of the mesh
  virtual dof_no_t nUnknowns();

protected:

  //! initializes the vectors and stiffness matrix with size
  virtual void createPetscObjects() = 0;

  const DihuContext context_;     ///< the context object with python config of the class that uses this data object
  std::shared_ptr<BasisOnMeshType> mesh_; ///< the mesh on which the data in this object is defined

  std::shared_ptr<Partition::RankSubset> rankSubset_;  ///< a subset of MPI ranks that will operate on the data of this object
  
  bool initialized_ = false;
};

} // namespace

#include "data_management/data.tpp"

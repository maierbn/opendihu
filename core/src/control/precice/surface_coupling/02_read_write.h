#pragma once

#include <Python.h>  // has to be the first included header

#include "control/precice/surface_coupling/01_initialize.h"

namespace Control
{

/** Generic Precice adapter, can be configured to either prescribe Neumann or Dirichlet boundary conditions.
 */
template<typename NestedSolver>
class PreciceAdapterReadWrite :
  public PreciceAdapterInitialize<NestedSolver>
{
public:

  //! constructor
  using PreciceAdapterInitialize<NestedSolver>::PreciceAdapterInitialize;

#ifdef HAVE_PRECICE
  //! read the incoming data from precice and set their values in the solver
  void preciceReadData();

  //! send the specified data over precice
  void preciceWriteData();
#endif

protected:

#ifdef HAVE_PRECICE
  //! set the data from preciceData as Neumann boundary conditions
  void setNeumannBoundaryConditions(typename PreciceAdapterInitialize<NestedSolver>::PreciceData &preciceData);

  //! set the data from preciceData as Dirichlet boundary conditions
  void setDirichletBoundaryConditions(typename PreciceAdapterInitialize<NestedSolver>::PreciceData &preciceData);
#endif

  std::vector<double> displacementValues_;   //< read value buffer for the displacement values
  std::vector<double> velocityValues_;      //< read value buffer for the velocity values
  std::vector<double> tractionValues_;        //< read value buffer for the traction values
  std::vector<Vec3> displacementVectors_;           //< send value buffer for the displacement values
  std::vector<Vec3> velocityVectors_;              //< send value buffer for the velocity values
  std::vector<Vec3> tractionVectors_;                //< send value buffer for the traction values

};

}  // namespace

#include "control/precice/surface_coupling/02_read_write.tpp"

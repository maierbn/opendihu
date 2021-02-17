#pragma once

#include <Python.h>  // has to be the first included header

#include "control/precice/volume_coupling/00_initialize.h"

namespace Control
{

/** Generic Precice adapter, can be configured to either prescribe Neumann or Dirichlet boundary conditions.
 */
template<typename NestedSolver>
class PreciceAdapterVolumeCouplingReadWrite :
  public PreciceAdapterVolumeCouplingInitialize<NestedSolver>
{
public:

  //! constructor
  using PreciceAdapterVolumeCouplingInitialize<NestedSolver>::PreciceAdapterVolumeCouplingInitialize;

#ifdef HAVE_PRECICE
  //! read the incoming data from precice and set their values in the solver
  void preciceReadData();

  //! send the specified data over precice
  void preciceWriteData();
#endif

protected:

  std::vector<double> scalarValues_;         //< temporary buffer for any scalar values
  std::vector<double> scalarValuesOfMesh_;   //< second buffer for any scalar values
  std::vector<Vec3> geometryValues_;         //< temporary buffer for the geometry values

};

}  // namespace

#include "control/precice/volume_coupling/01_read_write.tpp"

#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "control/precice/volume_coupling/01_read_write.h"

namespace Control
{

/** Generic Precice adapter, can be configured to either prescribe Neumann or Dirichlet boundary conditions.
 */
template<typename NestedSolver>
class PreciceAdapterVolumeCoupling :
  public PreciceAdapterVolumeCouplingReadWrite<NestedSolver>
{
public:

  //! define the type of the data object
  typedef typename NestedSolver::Data Data;

  //! constructor
  using PreciceAdapterVolumeCouplingReadWrite<NestedSolver>::PreciceAdapterVolumeCouplingReadWrite;

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary ("uninitialize")
  void reset();

  //! return the data object, with the call to this method the output writers get the data to create their output files
  Data &data();

};

}  // namespace

#include "control/precice/volume_coupling/precice_adapter_volume_coupling.tpp"

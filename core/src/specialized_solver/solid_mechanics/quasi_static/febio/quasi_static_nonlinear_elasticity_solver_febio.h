#pragma once

#include <Python.h>  // has to be the first included header

#include "data_management/specialized_solver/quasi_static_nonlinear_elasticity_febio.h"
#include "specialized_solver/solid_mechanics/quasi_static/febio/nonlinear_elasticity_solver_febio.h"
#include "output_connector_data_transfer/output_connector_data.h"
#include "output_writer/manager.h"

namespace TimeSteppingScheme
{

/** A specialized solver for 3D linear elasticity, as quasi-static timestepping scheme (a new static solution every timestep).
 *
 *  This class uses febio2.
  */
class QuasiStaticNonlinearElasticitySolverFebio :
  public NonlinearElasticitySolverFebio
{
public:

  //! constructor
  QuasiStaticNonlinearElasticitySolverFebio(DihuContext context);

protected:

  //! create the febio_input.feb file which contains the problem for febio to solve
  virtual void createFebioInputFile() override;

  double force_;                    //< factor of force that is applied in axial direction of the muscle
};

}  // namespace


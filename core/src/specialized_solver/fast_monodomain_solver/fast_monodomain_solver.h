#pragma once

#include <Python.h>  // has to be the first included header

#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

//#include "hodgkin_huxley.h"
//#include "tk2014.h"

/** \brief An efficient solver for the monodomain equation, as used in the fibers_emg example.
 *
 *  This is like a wrapper to an existing solver structure of Strang with Heun (0D problem) and Implicit Euler (1D problem).
 *  It solves the same equation as the wrapped solver would, but more efficiently.
 *  The possibly distributed fiber data is communicated to one rank and computed there, such that every fiber is always computed on a single rank.
 *  After all timesteps are done, the data is communicated back, such that from outside this solver appears equal to the corresponding nested solver structure.
 *  Internally, efficient vectorized code is produced using Vc and the CellmlSourceCodeGenerator.
 *
 *  Stimulation of fibers is done using a fiberDistributionFile and a firingTimesFile.
  */

/** This is the general class, relevant are the specializations.
  */
template<typename FibersEMG>
class FastMonodomainSolver
{
  static_assert(sizeof(FibersEMG) == 0, "This specialization of FastMonodomainSolver is not implemented. Currently only Hodgkin-Huxley and new_slow_TK are implemented!");
};

/** The main implemented class.
 */
template<int nStates, int nIntermediates>
class FastMonodomainSolver
<
  Control::MultipleInstances<                       // fibers
    OperatorSplitting::Strang<
      Control::MultipleInstances<
        TimeSteppingScheme::Heun<                   // fiber reaction term
          CellmlAdapter<
            nStates,nIntermediates,
            FunctionSpace::FunctionSpace<
              Mesh::StructuredDeformableOfDimension<1>,
              BasisFunction::LagrangeOfOrder<1>
            >
          >
        >
      >,
      Control::MultipleInstances<
        TimeSteppingScheme::ImplicitEuler<          // fiber diffusion, note that implicit euler gives lower error in this case than crank nicolson
          SpatialDiscretization::FiniteElementMethod<
            Mesh::StructuredDeformableOfDimension<1>,
            BasisFunction::LagrangeOfOrder<1>,
            Quadrature::Gauss<2>,
            Equation::Dynamic::IsotropicDiffusion
          >
        >
      >
    >
  >
> : public FastMonodomainSolverBase<nStates,nIntermediates>
{
public:
  using FastMonodomainSolverBase<nStates,nIntermediates>::FastMonodomainSolverBase;
};

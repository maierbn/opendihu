# Surface Coupling Examples

This folder contains examples were muscles and tendons are coupled using preCICE, a coupling library for partitioned multi-physics simulations. 
PreCICE does not provide an official adapter for OpenDiHu (yet). Benjamin Maier developed an opendihu precice adapter which is included in the opendihu source code. The existing opendihu consists of two independent adapters; one for surface coupling and another one for volume coupling.  

The examples in this folder use the built-in opendihu precice adapter for surface coupling.

### About the muscle participant:

**A multi-scale problem**
The muscle is a multi-scale problem with three main phenomena being modelled:
* sub-cellular processes at the sarcomeres: modelled as a 0D process and with `dt_0D = 0.5e-3` [ms]                      
* potential propagation at the fibers: modelled as a 1D process and with `dt_1D = 1e-3` [ms] 
* solid-mechanics (muscle contraction): modelled as a 3D process and with `dt_elasticity = 1e-1` [ms] 

**The structure of the muscle solver**
The following opendihu solver is used:
``` 
Control::Coupling<
  MonodomainSolver,
  MuscleContractionSolver<
      Mesh::StructuredDeformableOfDimension<3>
  >
>
```

The muscle solver consists of the `MonodomainSolver`, which models the 1d voltage propagation in the muscle fibers, and the `MuscleContractionSolver`, which models the 3d solid mechanics. 

> **Note**
> **Connected slots between  `MonodomainSolver` and `MuscleContractionSolver`**
> The `MuscleContractionSolver` receives the following slots:
> `slotNames":                    ["m1lda", "m1ldot", "m1g_in", "m1T", "m1ux", "m1uy", "m1uz"]`
> to-do

The `MonodomainSolver` is defined as follows:
``` 
using MonodomainSolver =
  FastMonodomainSolver<                               
    Control::MultipleInstances<                       // one instance per fiber 
      OperatorSplitting::Strang<
        Control::MultipleInstances<
          TimeSteppingScheme::Heun<                   // fiber reaction term
            CellmlAdapter<
              N_STATES, N_ALGEBRAICS,                 // depends on the cellml model
              FunctionSpace::FunctionSpace<
                Mesh::StructuredDeformableOfDimension<1>,
                BasisFunction::LagrangeOfOrder<1>
              >
            >
          >
        >,
        Control::MultipleInstances<
          FiberDiffusionSolver<                       // fiber diffusion term 
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
  >;
```

The idea behind the `MonodomainSolver` is to apply strang splitting to each fiber. The components of the strang splitting are the fiber reaction term and the fiber diffusion term. The fiber reaction term corresponds to the sub-cellular processes taking place at the sarcomeres. 

> **Note**
> The overall timestep of the strang splitting is `dt_splitting_0D1D`
> The optimal choice of timesteps for the strang splitting is `dt_1D = 2*dt_1D`

> **Note**
> **Connected slots between reaction term and diffusion term**
> to-do


### About the tendon participant:
to-do

### About the muscle-tendon example:
**Set-up**
- muscle: fixed on one end (z=0.0) and attached to the tendon in the other end (z= muscle initial length)
- tendon: fixed on one end (z= muscle initial length + tendon length) and attached to the muscle on the other end (z= muscle initial length)
- the spatial discretization of the muscle and tendon is chosen so that we have a matching mesh (we define same number of elements on the x and y direction)
- transfer of data: the muscle tendon sends the displacement and velocities to the traction and the muscle sends traction to the tendon.
- coupling scheme: serial-explicit

If you run the simulation you will see this looks quite good! 

> **Note**
> There are only some mismatches for the traction values at the edges of the coupling interface, but I believed this is due to the implementation of von Neumann boundary conditions in openDiHu and is not a worrying issue.

**Open Issues**
- Implicit coupling reaches the maximum number of iterations (100). Currently trying absolute convergence criterium instead of relative criterium in case this is due to an almost-zero denominator.
- As a consequence of the previous point we cannot have running simulations where the traction is send from the tendon to the muscle.
- If we replace the free end of the tendon by a traction bc this boundary condition is not reflected in the results. However, a single tendon with different traction boundary conditions at the ends was simulated without issues.


# About the muscle-tendon-muscle example:

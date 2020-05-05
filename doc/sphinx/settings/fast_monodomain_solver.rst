FastMonodomainSolver
======================

This is a very efficient implementation of multiple fibers where the monodomain equation is solved on each.
The fibers all have the same number of elements. The fibers are more or less parallel. 
The parallel partitioning is arbitrary, i.e. fibers can be subdivided to different processes and different fibers do not need to share processes.

All nodes of this geometry form a cuboid in index space :math:`(i,j,k)`, where the `k` index runs over the nodes of a fiber and :math:`(i,j)` specify the fiber in a 2D grid of fibers.
The partitioning is obtained by dividing this cuboid by axis-aligned plane cuts in all three dimensions.

The *FastMonodomainSolver* reuses nested solvers, as given below.

Usage
----------
The following shows the typical usage in the C++ source file. Note that only the number of states and intermediates, 4 and 9 in this example can be changed.

.. code-block:: c
  :linenos:

  FastMonodomainSolver<                        // a wrapper that improves performance of multidomain
    Control::MultipleInstances<                       // fibers
      OperatorSplitting::Strang<
        Control::MultipleInstances<
          TimeSteppingScheme::Heun<                   // fiber reaction term
            CellmlAdapter<
              4, 9,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
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
  >

The two template arguments of `CellmlAdapter`, the *number of states* and *number of intermediates* can be adjusted to fit the subcellular CellML model.
All other templates must appear exactly as given above.

The *FastMonodomainSolver* solves the same equations as the nested solver would (just as if lines 1 and 27 were not present). Also the python options are the same. This class extracts most of its settings from the settings of the nested solvers.
For exemplary complete python options, see the ``multiple_fibers_cubes_partitiong`` or ``fibers_emg`` examples.
Note that the definition of the fiber meshes and ranks is a bit more involved, study the example to learn how it works.

The improved performance is by roughly a factor of 10. The reason for this is that the 1D diffusion problem which is a tri-diagonal system gets solved serially and by a Thomas' algorithm of linear time complexity. For this to work, all values of a fiber are communicated to a single rank. This fiber is then solved completely on this one rank for all specified timesteps. Afterwards the values are communicated back to the original process. Consequently, the *FastMonodomainSolver* appear to surrounding solvers like its nested solvers with a cubes-like partitioning, but internally the fibers are not split across processors.

For the subcellular model, efficent code of the whole Heun scheme, using `"vc"`, is generated and executed. This is again faster than if `Heun` and `CellMLAdapter` are nested.




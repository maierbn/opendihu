.. _introduction:

Introduction
============

Functioning of the musculoskeletal system is driven by processes on different temporal and spatial scales. On the smallest, molecular scale, ion transfers through the muscle fiber membranes and molecular processes on the "subcellular" level drive the behaviour of electrophysiology and force generation. Action potentials propagate along muscle fibers and contraction takes place on the muscle scale.

The solution of multi-scale models require a combination of tailored numerical schemes. Time step widths and spatial mesh resolutions have to be chosen carefully to avoid instabilities and allow a computation in feasible run times. Different meshes with different dimensionalities have to be defined and interpolation between them is required. The simulation has to run on multiple processes, to efficiently exploit today's hardware and reduce the runtime to a minimum. Input and output has to be processed in useful and efficient data formats for different purposes, for convienient debugging as well as handling of large datasets with production runs. Configuration of numerical, model and algorithmic properties has to be organised in a handy way.

OpenDiHu fulfills the previously mentioned requirements that are imposed from a users perspective. From the developers perspective, a modularized code structure in ensured by heavily using C++ templates to abstract concepts like meshes, function spaces and solvers.

Overview
------------

OpenDiHu is a software framework that solves static and dynamic multi-physics problems, spatially discretized in 1D, 2D and 3D by the finite element method.
Our core design goals are **usability**, **performance** and **extensibility**.

Good **usability** refers to the process of configuring a given simulation with the python interface, running variants of the simulation and producing output in helpful file formats. Because of the contained python interpreter reconfiguration is possible at runtime. 

The **performance** goal is satisfied within the C++ core implementation. The data structures avoid expensive data copy and allow for vectorization. All structured grid functionality is designed for parallel execution. We closely build on the parallely efficient `PETSc <https://www.mcs.anl.gov/petsc/>`_ library without overhead. The code was successfully run on 27,000 cores on the supercomputer "Hazel Hen" in Stuttgart.

The framework is **extensible** to future models. It provides infrastructure to store and manipulate scalar and vector fields, handle input and output, assemble system matrices and interface various PETSc solvers. For example, there is support for `CellML <https://www.cellml.org/>`_, a format for interchanging systems of ordinary equations. 
The framework and its applications are constantly extended. However, there are stable releases.

Read `my thesis <https://arxiv.org/abs/2107.07104>`__ or `one of the papers <https://github.com/maierbn/opendihu#literature--how-to-cite>`__ for a better introductory text ;-)

.. Documentation
.. ------------------

.. Theory documentation can be found in the `doc/derivations/doc.pdf` document. 
.. Some developer hints for the core code can be found in `doc/documentation.rst`. 
.. Generally the C++ code contains a lot of comments and often it is useful to look directly into the code, to find out how something works. 


.. toctree::
   :maxdepth: 1
   :caption: Next pages to read
   
   user/features
   user/installation
   user/framework_structure
   user/getting_started
   user/examples
   

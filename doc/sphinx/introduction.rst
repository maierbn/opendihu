.. _introduction:

Introduction
============

Biophysical processes can happen on several scales. In the case of a skeletal muscle, there are ion transfers on a molecular scale, happening at points on the muscle fiber membranes. Action potentials propagate along whole muscle fibers and contraction takes place on the muscle scale. Such multi-scale processes require a combination of tailored numerical solution schemes. Time step widths and spatial mesh resolutions have to be chosen carefully to avoid instabilities and allow a computation in feasible run times. Different meshes with different dimensionalities have to be defined and interpolation between those is required. The simulation has to run on multiple processes, to efficiently exploit today's hardware and reduce the runtime to a minimum. Input and output has to be processed in useful and efficient data formats for different purposes, for convienient debugging as well as handling of large datasets with production runs. Configuration of numerical, model and algorithmic properties has to be organised in a useful way.

Opendihu fulfills the previously mentioned requirements, which are imposed from a users perspective. From the developers perspective, a modularized code structure in ensured by heavily using C++ templates to abstract concepts like meshes, function spaces and solvers.

Overview
------------

Opendihu is a software framework that solves static and dynamic multi-physics problems, spatially discretized in 1D, 2D and 3D by the finite element method.
Our core design goals are **usability**, **performance** and **extensibility**.

Good **usability** refers to the process of configuring a given simulation with the python interface, running variants of the simulation and producing output in helpful file formats. Because of the contained python interpreter reconfiguration is possible at runtime. 

The **performance** goal is satisfied within the C++ core implementation. The data structures avoid expensive data copy and allow for vectorization. All structured grid functionality is designed for parallel execution. We closely build on the parallely efficient `PETSc <https://www.mcs.anl.gov/petsc/>`_ library without overhead. The code was successfully run on 27,000 cores on the supercomputer "Hazel Hen" in Stuttgart.

The framework is **extensible** to future models. It provides infrastructure to store and manipulate scalar and vector fields, handle input and output, assemble system matrices and interface various PETSc solvers. For example, there is support for `CellML <https://www.cellml.org/>`_, a format for interchanging systems of ordinary equations. 
The framework and its applications are constantly extended. However, there are stable releases.

.. Documentation
.. ------------------

.. Theory documentation can be found in the `doc/derivations/doc.pdf` document. 
.. Some developer hints for the core code can be found in `doc/documentation.rst`. 
.. Generally the C++ code contains a lot of comments and often it is useful to look directly into the code, to find out how something works. 


.. toctree::
   :maxdepth: 1
   :caption: Next pages to read
   
   user/features
   user/framework_structure
   user/installation
   user/getting_started
   

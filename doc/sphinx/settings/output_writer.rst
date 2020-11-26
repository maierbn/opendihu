OutputWriter
=============

OutputWriters are used to write the computed simulation data to files. The data is organised in field variables, which are scalar or vector fields over the computational domain.

Note: To output matrices and right hand side vectors, e.g. in MATLAB compatible format, the direct options of the :doc:`solver` class should be used.

Several different output writers can be specified in the ``"OutputWriter"`` list.
The output writers are then called one after another and produce output in different formats, depending on their properties.

What kind of data, i.e. which field variables will be written, is fixed by the C++ code and cannot be adjusted in the python settings. In general, it contains all interesting field variables that are used during computation.

However, different output data can be written from output writers at different locations in the class hierarchy.

The ``"OutputWriter"`` keyword can appear in the settings for the following classes: 

* :doc:`cellml_adapter`
* :doc:`multiple_instances`
* :doc:`finite_element_method`
* :doc:`timestepping_schemes_ode` (Heun, ExplicitEuler, etc.)
* specialized solvers like :doc:`static_bidomain_solver`, :doc:`quasi_static_linear_elasticity_solver`

The following is an example of specifying output writers of all possible formats. This would output the same data in all possible file formats and as python callback. Multiple output writers of the same format would also be allowed.

.. code-block:: python

  "OutputWriter" : [
      {"format": "Paraview",   "filename": "out/filename", "outputInterval": 1, "binary": False, "fixedFormat": False, "onlyNodalValues": True, "combineFiles": False},
      {"format": "PythonFile", "filename": "out/filename", "outputInterval": 1, "binary": False, "onlyNodalValues": True},
      {"format": "ExFile",     "filename": "out/filename", "outputInterval": 1, "sphereSize": "0.005*0.005*0.01"},
      {"format": "MegaMol",    "filename": "out/filename", "outputInterval": 1},
      {"format": "PythonCallback", "callback": callback,   "outputInterval": 1}
    ]

The formats are explained in more detail below.
There are common properties, like filename and outputInterval.

filename
----------

The file name of the output file to write. It can contain a path to the output directory (e.g. *out/*). If this directory does not yet exist, it will be created. Usually, for dynamic problems, files will be written for multiple timesteps. The filename will get a consecutive number of 7 digits appended (with leading zeros), then the MPI rank number of the process, then the format specific suffix. For example the *out* folder could contain the following files for the Paraview output writer:

.. code-block:: python

  filename_0000000.0.vtr
  filename_0000000.1.vtr
  filename_0000000.pvtr
  filename_0000001.0.vtr
  filename_0000001.1.vtr
  filename_0000001.pvtr
  filename_0000002.0.vtr
  filename_0000002.1.vtr
  filename_0000002.pvtr

In this case, the ``*.pvtr`` is the master file that references the two files ``*.0.vtr``, written by rank 0 and ``*.1.vtr`` written by rank 1. All these files can be loaded in `Paraview <https://www.paraview.org/>`_ at once.

outputInterval
----------------
The interval in which timesteps an actual file should be written. This is an integer value. The output writer is usually called after each time step. However, it only outputs a file after ``outputInterval`` timesteps have passed. This is useful to reduce the number of output files. 

A common practice to compute outputInterval is to write

.. code-block:: python

  "outputInterval": int(1./dt*output_timestep)
  
where ``dt`` is the timestep width of the scheme in seconds and ``output_timestep`` is the output interval in seconds. Then files are generated every ``output_timestep`` seconds of simulation time. Note the cast to ``int`` to get an integer value.

fileNumbering
---------------
*Default: "incremental"*

Defines how the output files should be numbered. With ``"incremental"`` the files get incremental number suffixes starting from 0. With ``"timeStepIndex"`` the file suffix corresponds to the time step index.  This means that the suffixes are not incremental if ``outputInterval`` does not equal 1. The index is counted on a per-integrator basis. That means, that each time a time step is performed with a specific integrator, the index for that integrater increases.

Paraview
------------
`Paraview <https://www.paraview.org/>`_ is a postprocessing tool that can efficiently handle large data and can also be executed in parallel. It supports file formats that can also be handled by the `Visualization Toolkit <https://vtk.org/>`_ (*VTK*). The output files can be ASCII-based or binary. Separate files for every process or combined files can be written and parsed by Paraview.

A file can always only contain one mesh (with the exception of using the ``combineFiles`` option), for different meshes like, e.g. 1D and 3D meshes, separated files have to be written.

Example (normally the properties of a single output writer are written in a single line):

.. code-block:: python

  "OutputWriter" : [
    {
      "format": "Paraview", 
      "filename": "out/filename", 
      "outputInterval": 1, 
      "binary": False, 
      "fixedFormat": False, 
      "onlyNodalValues": True, 
      "combineFiles": False
    },
  ]

binary
~~~~~~~
*Default: True*

Whether to produce binary data files. The output files are valid XML-documents, the headers and describing structure of the file is therefore always human-readable. However, the data block can either contain the floating point values as ASCII-text or be in a binary format.
Binary files are encoded in `Base64 <https://en.wikipedia.org/wiki/Base64>`, an encoding that is 4/3 the size of the actual raw data but only uses printable characters and, thus, gives valid XML documents.

fixedFormat
~~~~~~~~~~~~
*Default: True*

If ``binary: False`` is set, the values will be written as ASCII text to the file. If ``fixedFormat`` is set to ``True``, the values will be in a column-fixed format with a length of 16 characters and in scientific notation.


combineFiles
~~~~~~~~~~~~~
*Default: False*

For serial execution, normally every process writes their own files. However, for massively parallel runs, this gives a high load on the file system. In theses cases, one should enable this option. Then a single file will be written per timestep, using MPI collective file operations.

The collective files will also gather all 1D, 2D and 3D meshes, respectively. This means that one file containing all 1D meshes will be created, another one containing only 2D meshes and another one with 3D meshes, if there are any. This is useful in a scenario of numerous 1D muscle fibers. Without this option, a new file would be created for every muscle fiber, because it is a new mesh. With this option, all fibers are contained in a single file.

File suffixes
~~~~~~~~~~~~~~
Depending on the :doc:`mesh`, different file formats with different file endings are created.

* ``*.vtr`` StructuredRegularFixed mesh, serial execution
* ``*.pvtr`` StructuredRegularFixed mesh, master file for parallel execution, besides this file, every rank writes a ``*.vtr`` file.
* ``*.vts`` StructuredDeformable mesh, serial execution
* ``*.pvts`` StructuredDeformable mesh, master file for parallel execution, besides this file, every rank writes a ``*.vts`` file.
* ``*.vtu`` Unstructured mesh, serial execution
* ``*.pvtu`` Unstructured mesh, master file for parallel execution, besides this file, every rank writes a ``*.vtu`` file.
* ``*.vtp`` A file containing multiple 1D meshes, when option ``combineFiles`` is used

PythonFile
----------
A ``*.py`` file will be written, that can be parsed using python and the *pickle* package. Human-readable output files and binary files can be created, depending on the ``binary`` option. Both options are valid pickle formats.

The file contents can directly be loaded using *pickle* or with the ``py_reader`` utility, that load the data of potentially multiple input files, regardless whether they are binary or not, and also stitching together the pieces from different processes. An example code is given here:

.. code-block:: python

  import py_reader

  filenames = ["..."]

  # load all files
  data = py_reader.load_data(filenames)

The data can also be inspected using the ``catpy`` and ``plot`` utilities, that come with opendihu and are located in the ``scripts`` subdirectory.

Using 

.. code-block:: bash

  plot *.py

all files in the current folder will be plotted using *Matplotlib*. An animation will be saved in ``mp4`` format and shown in an interactive plot.

The utility

.. code-block:: bash

  catpy *.py

outputs the contents of the files to the console. This is useful for debugging and works again regardless if the files are binary or ASCII-based.

The command line arguments to these two utilities are simply all files that should be considered, possibly from multiple timesteps and from multiple processes. The plot script only handles 1D and 2D plots, but automatically detects the contents and adjusts the plot format accordingly.

PythonCallback
---------------
This output writer does not write any files by itself. Instead, it calls a python callback function with an object of what would be contained in the file when the ``PythonFile`` format would be used. This callback function can then do whatever the user wants and write the data in a custom format.

ExFile
-------
The EX file format is an ASCII-based file format for unstructured meshes that is used by `OpenCMISS <http://opencmiss.org>`_. EX files are only suited for small problem sizes. It is also output by `OpenCMISS Iron <http://opencmiss.org>`_ and can be visualized using `CMGUI <http://physiomeproject.org/software/opencmiss/cmgui/download>`_.
More details on this file format can be found `here in the opencmiss documentation <http://opencmiss.org/documentation/data_format/ex_file_format.html>`_.

A geometry description consists of an *\*.exelem* file that contains element adjacency information and an *\*.exnode* file with the actual node positions.
The output writer also outputs a *\*.com* file, which contains *Perl* commands for *CMGUI*. You can visualize the data in the Ex files by calling 

.. code-block:: bash

  cmgui file.com

This creates predefined visualizations in CMGUI that can be further adjusted.

The ``sphereSize`` option defines how spheres, used to visualize nodes, will be rendered. The format is ``x*y*z`` and the default is ``0.005*0.005*0.01``.

MegaMol
--------

`MegaMol <https://megamol.org/>`_ is a visualization framework started by `VISUS <https://www.visus.uni-stuttgart.de/en>`_ at University of Stuttgart. It efficiently visualizes spheres, also on supercomputing hardware.

The MegaMol output writer outputs files in the `Adaptable Input/Output System 2 (ADIOS2) <https://adios2.readthedocs.io/en/latest/>`_ format. MegaMol can directly read this format. If the file is written to ``/dev/shm/``, *In-Situ* visualization is performed that completely avoids the disc to generate visualization output.

Since the file format is binary packed and self-descriptive, it is also suited for long-term storage of the data or for large simulation output in general. However, it cannot be directly visualization with e.g. Paraview.

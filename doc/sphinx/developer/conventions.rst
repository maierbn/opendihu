
Conventions and Tips
========================

Style guide
-------------

- Write code in style similar to what is already there: camelCase, all in english, enough comments, comment all member variables and methods.
- Use less abbreviations, instead long names for variables and methods. (Except maybe in comments)
- The rest is common sense.

Commonly used terms
-----------------------

In this section some expressions are listed that occur frequently in variable names.

- *Element*:
    An element consists of has several *nodes*.
- *Node*: 
    A node has several *dofs*, for Lagrange ansatz functions there is one dof per node, for cubic Hermite there are 2^D dofs per node
    (e.g. 1D: 2 dofs per node for value and derivative)
- *Dof*: 
    "Degree of freedom", an element of the function space that is spanned by the basis functions on the mesh. 
    Dofs are located on the nodes, every node has at least 1 dof or more.
    To define a scalar field on the mesh one has to specify a value for each dof. For a vector field with nComponents components, (e.g. the geometryField with nComponents=3) there are
    `nComponents*nDofs` values needed.
- *Unknown*: 
    A component of a dof. If only scalar fields are considered, `nUnknowns = nDofs`, otherwise `nUnknows = nDofs*nComponents`
    Note that the total number of unknowns for which the system has to be solved can be a multiple of the number of dofs when there are multiple components. 
    (E.g. 3D finite elasticity with 3-component displacement and 1-component pressure has 4 components per dofs, i.e. `nUnkowns = nDofs*4`)
- *No*:
    Abbreviation for "number", in the sense of a counter value (i.e. not meaning "total sum" but as in numbering)
- *n*:
    E.g. nUnknowns, nDofs: "number", but this time meaning the total number/count of items

Vectors that contain multiple elements/nodes/dofs or methods that work on such vectors have to specify if the entities are in the local or global numbering, refering to parallel execution, where a process knowns only its own local subdomain.

Consistent naming of local and global quantities is given by the following schemes:

- `<x>LocalWithGhosts`: The items are numbered in the local numbering, all local items including the local ghosts are contained. The non-ghost locals are the first entries, the ghost follow at the end.
- `<x>LocalWithoutGhosts`: The items are numbered in the local numbering and do not contain ghosts.
- `<x>Global`: This means global numbering, which may not be consecutive, because only local information may be known.

Examplary method names that use these conventions are `nNodesLocalWithGhosts()` or `nDofsGlobal()`.

Numbering Schemes
----------------------

1. Global natural numbering:
     starting at 0 at the front left bottom corner, then continuing in x-direction, then y-direction, then z-direction (if in 3D)
     This numbering is used for global input of Dirichlet boundary conditions or right hand side values
2. Local numbering,
     per partition, starting from 0 to nDofsLocalWithoutGhosts-1, only numbering the non-ghost local dofs, in the same order as the global natural numbering but in the local partition.
     Then from nDofsLocalWithoutGhosts to nDofsLocalWithGhosts-1 the ghost dofs are numbered, again in the order given by the global natural numering.
     This numbering is used for accessing local vectors and matrices. One can access all values including ghosts by using the whole range or only access the non-ghost values by using the range [0,nDofsLocalWithoutGhosts-1]
3. Global Petsc numbering,
     starts at 0 on rank 0 and follows the local numbering for the non-ghost dofs, then it continues on rank 1 and follows the non-ghost dofs there and so on.
     This is the numbering that has to be used for accessing global Petsc vectors and matrices, however this is not needed because one can access these vectors and matrices through the local vectors and matrices.
     It is needed for the creation of the Vecs and Mats.

Unit tests
--------------

The unit tests are located in `testing/unit_testing/src`. They are automatically compiled after the library by `scons`.
There are tests for 1 rank, 2 ranks and 6 ranks. All of the respective tests are compiled as one executable, which runs all tests and fails if any test fails.

The executable can be run manually e.g. in  testing/unit_testing/build_debug/1_rank_tests`. 
To only run a single test, use the `--gtest_filter=<test name>` command. It also supports wild cards (*), e.g.
  
.. code-block:: bash

  ./1_rank_tests --gtest_filter=LaplaceTest.*Structured* -v
  
This runs the tests.
When working on unit tests, you can temporarily enable just the test you are working on in the file `testing/unit_testing/SConscript`.
There you can set the file names to be included in the executable and fully enable/disable the 1_rank/2_ranks/6_ranks tests. This reduces compile time.

Howto debug
---------------

Using GDB 
^^^^^^^^^^
  It is convient to define the following alias:

  .. code-block:: bash
  
    alias gdbrun='gdb -ex=run --args '
    
  then simply run with gdb:
  
  .. code-block:: bash
  
    gdbrun ./executable ../settings.py <further-arguments>
    
  After the program crashes you can inspect the stacktrace with the command `bt` inside gdb.


Debugging output
^^^^^^^^^^^^^^^^^^
  Logging is done using the `Easylogging++ <https://github.com/zuhd-org/easyloggingpp>`_ library. In the C++ code the statements are, e.g., `LOG(DEBUG) << "text";` and `VLOG(1) << "level1";`.
  The `LOG(DEBUG)` statements are only displayed if compiled to debug target. The `VLOG(1)`, `VLOG(2)`, etc. statements are also only available in debug executable and have to be enable additionally using the following command line arguments:
  
  - run with `-v` to enable all verbose output
  - run with `--v=1` or `--v=2` etc. to enable verbose output to a given level
  - run with `-vmodule=partitioned_petsc_vec_structured.tpp=2,01_mesh_partition_structured.tpp=1` to enable verbose output of level2 only in the file `partitioned_petsc_vec_structured.tpp` and verbose output level 1 only in file `01_mesh_partition_structured.tpp`
    Also wildcards (*) can be used, e.g. `-vmodule=*vec*=3,*mat*=5,*mesh_partition*=1`, then all files matching `*vec*`, `*mat*` or `*mesh_partition*` will get the specified output verbosity.
  
Debugging parallel programs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  - run program with mpirun and with `-pause` argument, example:
    
    .. code-block:: bash
  
      mpirun -n 2 ./2_ranks_tests -v --gtest_filter=LaplaceTest.Structured1DHermite -pause
      
    Then it will stop with the following message:
    
      0/2 INFO : Rank 0, PID 16614 is waiting for gdbResume=0 to become 1 
    
      gdb -p 16614
      
      select-frame 2
      set var gdbResume = 1
      info locals 
      continue
      1/2 INFO : Rank 1, PID 16615 is waiting for gdbResume=0 to become 1 
      
      gdb -p 16615
      
      select-frame 2
      set var gdbResume = 1
      info locals 
      continue
      
  - now in two separate shell windows, execute `gdb -p 16614` and `gdb -p 16615`. This attaches gdb to the two MPI processes. Inside gdb run the displayed commands `select-frame 2`, `set var gdbResume = 1`, etc.
    After `continue` in both attached shells the program will continue. When it crashes, use `bt` to inspect the location again.
    
Using Memcheck
^^^^^^^^^^^^^^^^
  
  For segmentation faults that cannot be debugged with gdb, you can use valgrind with memcheck:
  
  .. code-block:: bash
  
    valgrind --tool=memcheck ./executable
    
  There are a lot of "false positives" at the beginning while the python settings script is run. This is due to the python library overloading functions of memory management. Watch out for errors after these big outputs.
  
Working with parallel vectors
---------------------------------

The objects that represent parallel vectors in opendihu are of type `FieldVariable`.
A field variable is a vector with one entry for each dof in the `FunctionSpace`.
Each entry has a number of components, so the actual "size" of a field variable is `nDofsGlobal*nComponents`.

Examples for the use of components are the geometryField, which stores x,y,z values for each dofs and thus has `nComponents=3`.

To access the values of a field variable, there are methods like:

===========================  ==========================
  `getValue`, `getValues`    get the values for specified dof numbers and components
  `getElementValues`         get all the values that correspond to the dofs of an element
  `getValuesWithGhosts`      get all locally stored values, including ghost values
  `getValuesWithoutGhosts`   same, but without ghost values
===========================  ==========================

(Read the actual signatures in `field_variable/structured/03_field_variable_set_get_structured.h`)

To set values in the vector, there are similar methods

  setValue, setValues, zeroEntries, setValuesWithGhosts, setValuesWithoutGhosts

It is most efficient to get/set multiple values at once instead of calling getValue/setValue for every single update.
The described methods work fully with the local dof numbering and only modify local or ghost values.

To work with these vectors using Petsc there is the `valuesGlobal()` method that returns a global Petsc vector that can be used with e.g. `MatMult(matrix->globalValues(), fieldVariable->globalValues(), result);` etc.
Do not use Petsc routines to get and set values (`VecGetValues`, `VecSetValues`) with the obtained Petsc vectors! 
Rather use the described setValues/getValues methods, because they take care of the correct
indexing (local vs. global numbering) and perform sanity checks for indices (of course only when compiled for debug target).

There is also a more low-level vector class, `PartitionedPetscVec` which wraps the Petsc Vec's and handles creation of the Vec's and ghost values exchange. This class is used by FieldVariable internally and there should be no need to use it directly. However, there you can see how numbering/ghost exchange etc. is implemented.

One important thing when working with field variables, i.e. parallel vectors, is the correct ghost value manipulation.
Each rank only has local memory to store its non-ghost and ghost values.

Before you can read and write to locally stored ghost values, call 

.. code-block:: c

  fieldVariable->startGhostManipulation()

This fills the local ghost value buffer with the actual values from the corresponding other ranks and overwrites what was previously in the local ghost value buffer. After that you can read out the ghost values and also write to the local buffer. Calls with INSERT_VALUES and ADD_VALUES can be mixed without further consideration, because everything is only updated locally. For example you could do
  
.. code-block:: c

  fieldVariable->setValues(<vector of local dof nos for some non-ghosts and ghosts>, <some values>, INSERT_VALUES);
  fieldVariable->setValues(<some other dof nos>, <some other values>, ADD_VALUES);
  fieldVariable->getValues(<again some dof nos with possibly ghosts>, <output vector>);

After that for each ghost dof the ghost value on the rank where it is a ghost and the actual value on the rank where it is not a ghost need to be added. This is done by 

.. code-block:: c

  fieldVariable->finishGhostManipulation()

After that the two values are added and stored only on the rank where the dof is not the ghost. To also get the updated value to the rank where it is a ghost you need to call fieldVariable->startGhostManipulation() again. For every  startGhostManipulation there has to be a matching finishGhostManipulation later.

Note, the following is wrong:

.. code-block:: c

  fieldVariable->startGhostManipulation()
  // setValues which manipulates local ghost values
  fieldVariable->finishGhostManipulation()  // everything good up to here, now every rank has the correct local non-ghost values

  fieldVariable->startGhostManipulation()  // still okay, now every rank also has correct ghost values (#)
  // setValues which only manipulate non-ghost values
  fieldVariable->finishGhostManipulation()  // unexpected result, some local values (those that are ghosts on other ranks) will get ghost-buffer values added, that are still in the ghost buffers on an other rank. (#)

The correction for the example would be to remove (#) lines or set the ghost buffers to zero fieldVariable->zeroGhostBuffer() (but then the start/finishGhostManipulation calls would be useless anyway)

So if you want to read ghost values, call startGhostManipulation() beforehand, 
if you want to write all ghost values, wrap the setValues code with startGhostManipulation() and finishGhostManipulation(). 
If you want to write some ghost values, call startGhostManipulation(), save the ghost values you need (by fieldVariable->getValues()), zeroGhostBuffer(), finishGhostManipulation()

Using output data
-----------------------

The python output data in `*.py` files can be viewed by the script `catpy.py <files>` and plotted by `plot.py <files>`. 
The scripts are located in the `scripts` folder.
It is convenient to add this folder to PATH, e.g. in `~/.bashrc` with 

.. code-block:: bash
  
  export PATH=$PATH:/store/software/opendihu/scripts   # (adjust to your path)
  
There are also shortcuts `plot` and `catpy`. Example usage:

.. code-block:: bash

  catpy                   # without arguments, displays contents of all files in the current directory with `*.py` suffix.
  plot                    # without arguments, plots everything from the `*.py`files in the current directory.
  plot out*               # plot all out* files
  validate_parallel.py    # without arguments, checks if the content of all files with corresponding names matches, where some files are serial files like `out.py` and some are parallel files like `out.0.py`, `out.1.py` etc.
  validate_parallel.py out.py out.0.py out.1.py   # do the same but with explicit specification of which files to use.

Profiling on Hazel Hen
-------------------------


In order to use `pat_run` with GCC do the following:

- It needs `-finstrument-functions` in compile options for all code levels. I.e., core, example code and eventually in cellml runtime compiled code.
- Load the modules: `PrgEnv-gnu, perftools-base, perftools-preload`.
- Execute the aprun command: aprun (aprun opts like -n= ...) pat_run (pat_run options, 
  at least: -gmpi -r -m lite-events. for more, see man pat_run) example_name (example options)
  for example: ´aprun -n8 pat_run -m lite-events -gmpi -r ./shorten_implicit ../settings.py´

Work with large simulation output on servers
----------------------------------------------

Assume you have a long simulation on a compute server that produced a lot of output files. 
Now you want to download them to your laptop and visualize them, but the total filesize is too high.

- One solution is to use start the render server of paraview on the compute server and then connect from the local system with the Paraview client. 
  This requires a good internet connection, otherwise the interactivity is reduced. 
  
  Download and build paraview yourself, e.g., on neon. On your laptop, add a server launch configuration (Server Type `Client / Server`, Host `localhost`, Port `11116`, configure Startup Type `Command` and insert a command like the following (adjust to your path):
  
  .. code-block:: bash
  
    ssh -X maierbn@neon /home/maierbn/software/ParaView-5.6.0-RC1-Qt5-MPI-Linux-64bit/bin/pvserver --server-port=11116
    
  Then you just need to click on `Connect` in ParaView on your laptop. This will automatically run the server on neon. (You need to be in the network or VPN).
   
- Or you only select a subset of the files on the server to download them. Selecting every `nth` file can be done in bash:

  .. code-block:: bash
  
    tar czf ramp.tgz fibers_0000*{000..286..8}.vtp
    
  This selects every 8th file out of the files from `fibers_0000000.vtp` to `fibers_0000286.vtp` and puts them in a compressed archive, which can then be downloaded.


Building the  sphinx doc (what you are reading)
-------------------------------------------------

The following has to be installed beforehand.

.. code-block:: bash

  sudo pip install sphinx
  sudo pip install recommonmark
  sudo pip install sphinx_rtd_theme
  
To build the newest documentation, change into the opendihu/doc/sphinx directory and run

.. code-block:: bash

  make html

The local documentation can be viewed in a browser at

  file://<your-path>/opendihu/doc/sphinx/_build/html/index.html

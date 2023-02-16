
Getting started
=====================

* To get started, you'll find some examples in the `examples` directory. Also the system tests under `testing/system_testing/tests` might be useful to look at.

* To build an example, `cd` into a subdirectory under `examples`, e.g. `examples/laplace/laplace2d`. In this directory run ``scons``. 
  For this to work, you either need to install `scons` on your system (e.g. `sudo apt install scons` on ubuntu). Or use the given `scons` in the `dependencies` directory, ``../../dependencies/scons/scons.py``. See the `:doc:installation` page for details.

* To build the release target, simply use the alias ``sr`` (if defined, see :ref:`installation_aliases`) or use ``scons``, ``scons BUILD_TYPE=release`` or ``scons BUILD_TYPE=r``, to build the debug target, use ``sd`, ``scons BUILD_TYPE=debug`` or ``scons BUILD_TYPE=d``.

* To use multiple processes at once, add the ``-j`` option, e.g. ``scons BUILD_TYPE=r -j 4`` to build the release target with 4 processes. This is not needed for the ``sr`` shortcut as this does it automatically.

  .. code-block:: bash

    cd $OPENDIHU_HOME/examples/laplace/laplace2d
    mkorn && sr

* There will be executables created in the `build_debug` or `build_release` subdirectories. Change into one of these directories and run the program with a settings file as only argument: 

  .. code-block:: bash

    cd build_release
    ./laplace_regular settings_lagrange_quadratic.py

* Output files in this example (and likewise in the other examples) will be created under the `out` subdirectory. 
  If you look into `out`, you'll find two files: `laplace.py` and `laplace.vtr`.
 
  The `*.vtr` files can be visualized using Paraview. The `*.py` files can be visualized using the plot script in `opendihu/scripts`. 
  
  It is useful to add the `scripts` directory to the `PATH` environment variable, e.g. by adding the line `export PATH=$PATH:$OPENDIHU_HOME/scripts` to your `.bashrc` file as mentioned on the :doc:installation: page.
  Then you can run ``plot.py <*.py-output files>`` anywhere to visualize the output. Usually the shortcut ``plot`` instead of ``plot.py`` should also work. (Unless you have installed something different with the name ``plot``).
  
  An often used command is thus ``plot out/*``. If arguments are ommited, i.e. just ``plot``, it is the same as ``plot *.py``.
  In our example the command could be `
  
  .. code-block:: bash

    plot out/laplace.py
    
* The source files are located in the `src` subdirectory. The files to be compiled are specified in `SConscript`.
  There, you can, for example, comment out the examples that you don't want to compile every time.
* Now, change into `src` (i.e. into `.../examples/laplace/laplace2d/src`) and open `laplace_regular.cpp`. This is the main file of the 2D Laplace model. As can be seen, the equation `Δu = 0` is discretized by the Finite Element Method on a structured regular grid of dimension 2, basis functions are Lagrange of order 2, and the Quadrature scheme is Gauss quadrature with 3 gauss points per dimension. 
  
  * Now change the quadratic to linear Lagrange basis functions, by changing ``LagrangeOfOrder<2>`` to ``LagrangeOfOrder<1>``.
  * Change into the `build_debug` directory (``cd ../build_debug``). 
  * If you have set the aliases, you can recompile with ``sdd``. Otherwise go up one directory and run ``scons BUILD_TYPE=d -j 4``. 
  * Now, from the `build_debug` directory, run the new executable with 

  .. code-block:: bash

      ./laplace_regular ../settings_lagrange_linear.py
      
  * Plot the result with ``plot out/*``.

* Test the parallel execution and run the same program with the same settings file on two processes:

  .. code-block:: bash

    mpirun -n 2 ./laplace_regular ../settings_lagrange_linear.py

  If you now look into the out directory (``ls -l out``), you'll see that two new files, `laplace.0.py` and `laplace.1.py`, were created by the two processes. The file `laplace.py` is still the old one from the single process.

  Now plot the new files, either ``plot out/laplace.0.py out/laplace.1.py`` or shorter ``plot out/laplace.*.py``. The result is the same.

  Check that the results from the serial and parallel are actually the same using the following helper script:

  .. code-block:: bash

      validate_parallel.py out/*
      
* The created python output files are human-readable (because `"binary":False` is set in the settings file). You can open them in an editor and see what they contain. There is also the `catpy`  script for formatted printing on the console:

  .. code-block:: bash

    catpy out/laplace.0.py
    
* With the current settings, also the Paraview files are human-readable. You can also open e.g. `out/laplace.vtr` in an editor. Also try loading the `.pvtr` file in Paraview. 
  For big files it is better to produce binary files.
  
  In the settings file `settings_lagrange_linear.py` change `"binary":False` to `"binary":True` in the output writers. Now if you run the program again you'll get binary files that can't be read in a text editor. However, the ``plot``, ``validate_parallel`` and ``catpy`` utilities still work. 
* If you know `cmgui <http://physiomeproject.org/software/opencmiss/cmgui/download>`_, the visualization tool of `OpenCMISS <http://opencmiss.org/>`_ Zinc, you can also generate `exnode` and `exelem` output files for cmgui. Add the line

  .. code-block:: python

      {"format": "Exfile", "filename": "out/laplace"},
    
  to the `"OutputWriter"` list in file `settings_lagrange_linear.py` (line 31). (More details at :doc:`/settings/output_writer`.)
  After running the program again, you get the output files `laplace.exelem`, `laplace.exnode` and `laplace.com` in the out directory. The `.com` file is a convienient perl script that sets up the visualization in cmgui (OpenCMISS Iron won't generate this for you.). Change into the `out` directory and simply run `cmgui laplace.com`. In the Scene Editor click on `/` and then the `surface` item. Under `data`, select `solution` as the field variable that will be shown in color. Now you can tilt the view in the Graphics window to see the solution.
    
* Now you know the basics, how to run a simulation program. Next, you can try to change parameters in the settings file, like number of elements (variables `m` and `n`), the `physicalExtent` or try to understand, how the Dirichlet boundary conditions were specified. 
  Note, that because this example uses a `Mesh::StructuredRegularFixedOfDimension<2>` mesh (in the `cpp` source file), we can only have elements with quadratic shape, i.e. `physicalExtent` and `nElements` have to match. You can look into the `laplace_structured.cpp` example file, which uses a structured mesh, that can have different mesh width in `x` and `y` direction or even arbitrary node positions.
* The settings files use python syntax and are actually Python scripts. 
  This means you can execute any Python code there, for example load your own custom geometry or input data files and set the options appropriately. 
  The general documentation of the options is given on the :doc:`/settings` pages, 
  but some classes are not yet documented and their settings can only be known from the examples (or the C++ core code).
  So if you need to know how to specify certain options, look for an example files, that does it, or ask a developer.
* To execute some of the more advanced electrophysiology examples, you'll need special input files like a muscle geometry. 
  These are too large to have in git. `Download the files <https://zenodo.org/record/4705982>`_ and put them in the ``examples/electrophysiology/input`` directory.
* If you now continue to use opendihu, you can read the :doc:`/settings` pages for reference. 
  If anything is unclear do not hesitate to ask. If you have improvements concerning the formulations on this website or can contribute to writing the documentation, come in contact!

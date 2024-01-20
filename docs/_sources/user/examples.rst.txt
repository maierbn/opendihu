Examples
===========

The actual simulation programs or scenarios that make use of opendihu as software framework are called examples. They are stored in the `examples` folder, organized in subdirectories.

Directory structure
-------------------------

An example consists of the following files, as illustrated using the `electrophysiology/cellml` example.

.. code-block:: bash

  $ tree -L 2 cellml/
  cellml/
  ├── build_debug
  │   └── cellml
  ├── build_release
  │   └── cellml
  ├── SConscript
  ├── SConstruct
  ├── settings.py
  └── src
      └── cellml.cpp

The following explanation starts from the bottom. 

* The `src` subdirectory contains the source file of the example. This is the short C++ file, that instantiates the solver and runs the simulation. In this case, `src/cellml.cpp` is only 22 lines, other examples are slightly longer.
* `settings.py` is a python script that specifies the settings to the simulation. The script will be parsed by the simulation program, the settings have to be present in the `config` variable, a Python dict.
* `SConstruct` is the file that tells `scons` how to build everything, similar to a `Makefile` for `make`. Every example has a copy of this file. Note the line
  
  .. code-block:: python
  
    # get the directory where opendihu is installed (the top level directory of opendihu)
    opendihu_home = os.environ.get('OPENDIHU_HOME') or "../../.."

  This determines the home directory of the opendihu installation. It tries to read the path from the environment variable `OPENDIHU_HOME`. If this fails, because the variable is not set, it uses the relative path `../../..`. This has the following consequences. First, if you want to have your own examples stored somewhere else and their `SConstruct` file looks the same, because you copied it from an other example, you should set the environment variable `OPENDIHU_HOME` such that the base directory of opendihu can be found. But a line similar to the following to your `~/.bashrc`.
  
  .. code-block:: bash
  
    export OPENDIHU_HOME=/your-absolute-path-to/opendihu

  Second, if you create a new example under the examples folder, make sure that the relative path match the depth of the example directory, e.g. `../../..` or `../..`.
* The `SConscript` file contains the executable names and the corresponding source files. You can define multiple executables.
  The following is a schematic example for two executables with names `executable-name` and `other-executable` and different source files:

  .. code-block:: python
  
    env.Program(target = "executable-name", source = ["src/source.cpp"])
    env.Program(target = "other-executable", source = ["src/other_source.cpp"])

* The folders `build_debug` and `build_release` will contain the compiled executable. They can be deleted, because after compilation they will be created again.

To compile an example, change into its base directory (where `SConstruct` is located) and execute

.. code-block:: bash

  sd
  
to build the debug target or

.. code-block:: bash

  sr
  
to build the release target. This will create the `build_debug` or `build_release` subdirectories if they didn't exist yet. (Note, you need the aliases `sd` and `sr` as described in :ref:`installation_aliases`).

Create your own example
-------------------------

Copy any existing example and delete everything except `SConscript`, `SConstruct`, `settings.py` and the `src` folder. Adjust these files depending on your need. Notice the comments above concerning finding of the opendihu home directory.
  
Input files
-------------

Sometimes, additional input files are required. If they are big, they are not included in the git repository, but hosted somewhere else. This is the case for all electrophysiology examples. These have an `opendihu/examples/electrophysiology/input` directory with binary data. There is a readme file that contains a link where to download these files.

For a description of existing example, see :doc:`existing_examples`.

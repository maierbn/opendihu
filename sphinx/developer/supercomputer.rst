Build on supercomputer
=========================

On Hawk, we have the problem that there is no direct internet access from the login nodes and compute nodes.
This makes the installation more difficult.

Git tunnel
------------

Acess to git repositiories, however, is possible by setting up the ssh tunnel:

ssh -N -D 1080 localhost

add the following to your bashrc (replace `iccbnmai` with your username):

.. code-block:: bash
  
  alias hawk='echo "start ssh-tunnel with ssh -N -D 1080 localhost" && ssh -R 7480:localhost:1080 icbbnmai@hawk-login04.hww.hlrs.de'

Then, on the computer that was specified during admission of the account, execute:

.. code-block:: bash
   
   ssh -N -D 1080 localhost


And in another shell:
   
.. code-block:: bash

   hawk

On hawk, add the following to your `~/.gitconfig`:

.. code-block:: bash

  [http]
    proxy = socks5://localhost:7480

(Make sure the port number matches, in this case `7480`, the same port can be used only once, so maybe use a different number if this port is already occupied.)

Modules
---------
To install opendihu on Hawk, load at least the following modules:

.. code-block:: bash

  module load adios2/2.5.0
  module load cmake
  module load python
  module load mkl
  module load zlib  # (needed to build python)


Installation
----------------
After cloning opendihu, you have to download all required packages on a different computer and transfer them to Hawk.

Then, you can in principle run ``scons no_tests=TRUE -j 64`` on a login node and the dependencies will be installed. 
If it says something like `Downloading...`, then interrupt the process and see what is missing.

For some packages, you have to do the installation manually. 
If you're done run ``touch scons_build_success`` in the ``dependencies/package/src/package`` directory, such that scons knows that the package was installed.

In the following, we give some hints how to install some packages. However, you have to figure out the details yourself, adjust the paths, etc.

* Hints to manually install PETSc:

  Download the dependencies for PETSc, e.g., to `dependencies/petsc_downloads`:
  
  .. code-block:: bash
  
    total 41M
    -rw------- 1 icbbnmai cbm44102 5,5M Apr 18  2020 93baaa8c9.tar.gz
    drwxr-xr-x 4 icbbnmai cbm44102 4,0K Apr 19  2020 bison
    drwxr-xr-x 4 icbbnmai cbm44102 4,0K Apr 19  2020 flex
    -rw------- 1 icbbnmai cbm44102 5,8M Apr 18  2020 scotch-v6.0.8.tar.gz
    -rw------- 1 icbbnmai cbm44102 5,8M Apr 18  2020 scotch-v6.0.8.tar.gz.1
    -rw------- 1 icbbnmai cbm44102 8,0M Apr 18  2020 sundials-2.5.0p1.tar.gz
    -rw------- 1 icbbnmai cbm44102 4,6M Apr 18  2020 v2.0.2-p2.tar.gz
    -rw------- 1 icbbnmai cbm44102 4,6M Sep 15  2020 v2.1.0-p1.tar.gz
    -rw------- 1 icbbnmai cbm44102 3,2M Apr 18  2020 v3.4.2-p2.tar.gz
    -rw------- 1 icbbnmai cbm44102 143K Apr 18  2020 v4.0.3-p5.tar.gz
    -rw------- 1 icbbnmai cbm44102 281K Apr 18  2020 v5.1.0-p7.tar.gz
    -rw------- 1 icbbnmai cbm44102 2,9M Apr 18  2020 v5.2.1-p2.tar.gz

  Building somehow only works on the login node, not on the compute node.

  .. code-block:: bash

    export PATH=$PATH:/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2//opendihu-hawk-gnu/dependencies/petsc_downloads/flex/install/bin
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2//opendihu-hawk-gnu/dependencies/petsc_downloads/flex/install/lib
    export C_INCLUDE_PATH=$C_INCLUDE_PATH:/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2//opendihu-hawk-gnu/dependencies/petsc_downloads/flex/install/include
    export PATH=$PATH:/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2//opendihu-hawk-gnu/dependencies/petsc_downloads/bison/install/bin
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2//opendihu-hawk-gnu/dependencies/petsc_downloads/bison/install/lib
    export C_INCLUDE_PATH=$C_INCLUDE_PATH:/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2//opendihu-hawk-gnu/dependencies/petsc_downloads/bison/install/include

    module load parmetis/4.0.3-int32-shared mkl/19.1.0 metis/5.1.0-int32-shared scalapack/2.1.0-shared scotch/6.0.9-int32-shared mumps/5.2.1-int32-shared sundials/5.1.0-int32

    ./configure CC=$CC CXX=$CXX FC=$FC F77=$F77 F90=$F90 --prefix=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2//opendihu-hawk-gnu/dependencies/petsc/install  --with-shared-libraries=1 --download--fblaslapack=1 --with-packages-download-dir=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2//opendihu-hawk-gnu/dependencies/petsc_downloads  --download-hypre \
     --with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native'

    # to add mumps: --download-mumps --download-scalapack 
    # but does not compile

    make PETSC_DIR=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2//opendihu-hawk-gnu/dependencies/petsc/src/petsc-3.13.1 PETSC_ARCH=arch-linux-c-opt all
    make PETSC_DIR=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2//opendihu-hawk-gnu/dependencies/petsc/src/petsc-3.13.1 PETSC_ARCH=arch-linux-c-opt install

* Hints to manually install python:

  Run the normal configure with the correct ``--prefix <opendihu>/dependencies/python/install/bin`` option. Then manually install the python packages `numpy, scipy and matplotlib` as follows:
  Download from pip:
  
  .. code-block:: bash
  
    argparse-1.4.0-py2.py3-none-any.whl
    cycler-0.10.0-py2.py3-none-any.whl
    geomdl-5.3.0-py2.py3-none-any.whl
    kiwisolver-1.3.1-cp39-cp39-manylinux1_x86_64.whl
    matplotlib-3.4.2-cp39-cp39-manylinux1_x86_64.whl
    numpy-1.20.3-cp39-cp39-manylinux_2_12_x86_64.manylinux2010_x86_64.whl
    Pillow-8.2.0-cp39-cp39-manylinux1_x86_64.whl
    pyparsing-2.4.7-py2.py3-none-any.whl
    python_dateutil-2.8.1-py2.py3-none-any.whl
    pytz-2020.4-py2.py3-none-any.whl
    scipy-1.6.3-cp39-cp39-manylinux1_x86_64.whl
    scons_build_success
    setuptools-50.3.2-py3-none-any.whl
    six-1.15.0-py2.py3-none-any.whl

  Then run ``<opendihu>/dependencies/python/install/bin/python -m pip install *.whl``.
  
* Hints to manually install preCICE:

  .. code-block:: bash

    cmake \
      -DCMAKE_INSTALL_PREFIX=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2/opendihu-hawk-gnu/dependencies/precice/install \
      -DCMAKE_BUILD_TYPE=RELEASE \
      -DPYTHON_EXECUTABLE=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2/opendihu-hawk-gnu/dependencies/python/install/bin/python3 \
      -DPETSC_DIR=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2/opendihu-hawk-gnu/dependencies/petsc/install \
      -DPETSC_EXECUTABLE_RUNS=TRUE \
      -DPRECICE_ENABLE_FORTRAN=OFF \
      -DPETSC_COMPILER=mpic++ \
      -DEigen3_ROOT=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2/opendihu-hawk-gnu/dependencies/precice/src/precice-2.1.0/eigen-3.3.8 \
      -DLIBXML2_INCLUDE_DIR=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2/opendihu-hawk-gnu/dependencies/precice/install/libxml2/include/libxml2 \
      -DLIBXML2_LIBRARY=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2/opendihu-hawk-gnu/dependencies/precice/install/libxml2/lib/libxml2.so \
      -DBOOST_ROOT=/lustre/hpe/ws10/ws10.2/ws/icbbnmai-opendihu2/opendihu-hawk-gnu/dependencies/precice/install \
      ..
    make -j 64
    make install

* After having installed the dependencies, run ``make clean`` to clear the scons cache and ``scons no_tests=TRUE -j 64`` to build.
  Set ``<PACKAGE>_DOWNLOAD = False`` in `user-variables.scons.py` to disable the automatic building efforts of scons where have installed it manually.
  
  
  

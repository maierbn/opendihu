Validation
=======================

On this page, we show some numeric experiments to validate the implemented solvers.

1D Poisson 
--------------

We solve the 1D Poisson equation with linear, quadratic and cubic Hermite elements and conduct a convergence study.
The resulting linear system of equations is solved directly using LU decomposition. Note that Hermite ansatz is only experimental, 
while linear and quadratic ansatz functions can safely be used.

Problem definition
^^^^^^^^^^^^^^^^^^^^^^

The 1D Poisson problem on a finite interval is given by

.. math::

  -\dfrac{\partial}{\partial x}u(x) + f(x) = 0 \quad \text{on }\Omega = [0, b], \quad b = 3.

We use the following right hand side.

.. math::

    f: [0, b] \to \mathbb{R}, \quad f(x) := 1 - x^2

We impose Dirichlet boundary conditions at both ends of the interval.

.. math::

    u(0) = 1,\\
    u(3) = 2

Analytic solution
^^^^^^^^^^^^^^^^^^^^^^

The following function fulfills the formulated boundary value problem.

.. math::

    u(x) = -\dfrac1{12}x^4 + \dfrac1{2}x^2 + \dfrac{13}{12}x + 1

Numeric solution 
^^^^^^^^^^^^^^^^^^^^^

We discretize the 1D domain with linear, quadratic and cubic Hermite elements with a given numbers of elements, :math:`n`. 
The following plots show the result for :math:`n=6`.

.. _poisson_1d_linear6:
.. figure:: /user/validation/poisson_1d_linear6.png
  :width: 70%

  Numeric solution with 6 linear elements
  
.. _poisson_1d_quadratic6:
.. figure:: /user/validation/poisson_1d_quadratic6.png
  :width: 70%

  Numeric solution with 6 quadratic elements
  
.. _poisson_1d_hermite6:
.. figure:: /user/validation/poisson_1d_hermite6.png
  :width: 70%

  Numeric solution with 6 Hermite elements  

Convergence study
^^^^^^^^^^^^^^^^^^^^^
To assess the convergence behaviour, we increase :math:`n` from 10 to 100 and compute the :math:`\mathcal{L}_2`-error against the analytic solution.

.. _poisson_1d_convergence_study:
.. figure:: /user/validation/poisson_1d_convergence_study.png
  :width: 100%

  Result of a convergence study, L2 error against number of elements

The experimental order of converge :math:`\frac{log(\varepsilon_i)-log(\varepsilon_{i-1})}{log(n_i)-log(n_{i-1})}` is determined to be as follows.

.. list-table:: Experimental order of convergence for `poisson_1d_convergence_study`_.
   :widths: 25 25

   * - Linear elements
     - -2
   * - Quadratic elements
     - -4
   * - Hermite elements
     - -1

The linear and quadratic formulations exhibit the expected order of convergence of exactly second order for linear elements and at least third order for quadratic elements.
Note that the solution :math:`u(x)` has no :math:`x^3`-term, therefore the order of convergence is exactly 4 with the quadratic ansatz.

For the quadratic ansatz, the error does not continue to decrease further below :math:`2^{-11}` for :math:`n=600` elements and more because machine precision is reached in the solver and/or the computation of the L2 error.

The Hermite formulation in OpenDiHu is experimental. 
The observed convergence order lies around 1. Using Hermite ansatz functions required more effort in the problem specification, because we have to define also boundary conditions for the first derivative :math:`u'(0)`.
In OpenDiHu, using the Hermite ansatz requires specification of both the right hand side function :math:`f` and its derivative :math:`f'`.

How to reproduce
^^^^^^^^^^^^^^^^^^^^^

* Build the example

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/poisson_1d
        mkorn && sr       # build

* Run the simulation for a given number of elements, e.g., 6 elements.

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/poisson_1d/build_release
        
        ./linear ../settings_poisson_1d.py 6
        ./quadratic ../settings_poisson_1d.py 6
        ./hermite ../settings_poisson_1d.py 6

* Plot the results (make sure that the ``$OPENDIHU_HOME/scripts`` directory is in the ``PYTHONPATH`` environment variable).

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/poisson_1d/build_release/out
        plot linear.py
        plot quadratic.py
        plot hermite.py

        # (instead you can also call `plot out/linear.py` etc. from one directory above)

* Compute the error for the current simulation results.

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/poisson_1d
        python3 compute_error.py

* Conduct the convergence study.
  
    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/poisson_1d/build_release
        python3 ./do_convergence_study.py


3D Poisson
---------------
The next validation scenario is a more complex Poisson problem in a 3D domain.


Problem definition
^^^^^^^^^^^^^^^^^^^^^^

The 1D Poisson problem is given by

.. math::

  -\Delta u(\textbf{x}) + f(\textbf{x}) = 0 \quad \text{on }\Omega = [0, 2] \times [0,3] \times [0,4].

We define the following right hand side.

.. math::

    f: \Omega \to \mathbb{R}, \quad f(x,y,z) := 2\,x^3\,z + 24\,x^2\,y^2\,z + 8\,x^2\,z^3 + 6\,x\,y^2\,z + 12\,x\,y\,z + 8\,y^2\,z^3 + 4\,y

We use the following Dirichlet boundary conditions :math:`u(\textbf{x}) = \bar{u}(\textbf{x})` on :math:`\textbf{x} \in \partial \Omega`

.. math::

    \bar{u}(x,y,z)\vert_{x=0} &= 4\,y\,(2\,y\,z^3 + 1) +1\\
    \bar{u}(x,y,z)\vert_{x=2} &= 8\,y^2\,z^3 + 108\,y^2\,z + 24\,y\,z + 4\,y + 32\,z^3 + 16\,z +1\\
    \bar{u}(x,y,z)\vert_{y=0} &= 2\,x^2\,z\,(x + 4\,z^2)+1\\
    \bar{u}(x,y,z)\vert_{y=3} &= 2\,x^3\,z + 8\,x^2\,z^3 + 216\,x^2\,z + 90\,x\,z + 72\,z^3 + 13\\
    \bar{u}(x,y,z)\vert_{z=0} &= 4\,y+1\\
    \bar{u}(x,y,z)\vert_{z=4} &= 8\,x^3 + 96\,x^2\,y^2 + 512\,x^2 + 24\,x\,y^2 + 48\,x\,y + 512\,y^2 + 4\,y+1
    
Analytic solution
^^^^^^^^^^^^^^^^^^^^^^

The following function fulfills the formulated boundary value problem.

.. math::

    u(x) = x^3\,y^2\,z + 4\,x^2\,y^2\,z^3 + 2\,x\,y^3\,z - y\,z^2 + 3\,x^2\,y + 1

Numeric solution 
^^^^^^^^^^^^^^^^^^^^^

    We discretize the domain :math:`\Omega = [0, 2] \times [0,3] \times [0,4]` by 
    
    * \(a\) :math:`\quad 2\,n\times 3\,n \times 4\,n\quad` linear and quadratic elements, and by
    * \(b\) :math:`\quad n\times n \times n\quad` linear and quadratic elements.

    In case (a), the elements have a cube shape while in (b) they are cuboid shaped with different element lengths along the coordinate axes.

    The resulting solution can be seen in the three following plots, for :math:`n=2` and :math:`z=1, z=2`, and :math:`z=3`.

.. |poisson_3d_z1| image:: /user/validation/poisson_3d_z1.png
    :width: 30%

.. |poisson_3d_z2| image:: /user/validation/poisson_3d_z2.png
    :width: 30%

.. |poisson_3d_z3| image:: /user/validation/poisson_3d_z3.png
    :width: 30%

|poisson_3d_z1| |poisson_3d_z2| |poisson_3d_z3|

The resulting :math:`\mathcal{L}_2`-errors to the analytic solution are given in :numref:`eoftable`. It can be seen that the error is at machine precision for every tested discretization.

.. _eoftable:
.. list-table:: Experimental order of convergence for `poisson_1d_convergence_study`_
   :header-rows: 1

   * - case
     - ansatz functions 
     - number of degrees of freedom
     - :math:`\mathcal{L}_2`-norm of error
   * - \(a\)
     - Linear
     - 60
     - 1.1e-14
   * - 
     - 
     - 315
     - 8.7e-14
   * - 
     - 
     - 910
     - 1.9e-13
   * - 
     - 
     - 1989
     - 3.7e-13
   * - \(a\)
     - Quadratic
     - 315
     - 1.7e-13
   * - 
     - 
     - 1989
     - 6.4e-13
   * - 
     - 
     - 6175
     - 1.3e-12
   * - 
     - 
     - 14025
     - 1.3e-11
   * - \(b\)
     - Linear
     - 27
     - 2.0e-14
   * - 
     - 
     - 125
     - 6.9e-14
   * - 
     - 
     - 512
     - 1.7e-13
   * - 
     - 
     - 1989
     - 3.7e-13
   * - \(b\)
     - Quadratic
     - 125
     - 9.5e-14
   * - 
     - 
     - 729
     - 3.6e-13
   * - 
     - 
     - 3375
     - 1.8e-12
   * - 
     - 
     - 12167
     - 8.1e-12

How to reproduce
^^^^^^^^^^^^^^^^^^^^^

* Build the example

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/poisson_3d
        mkorn && sr       # build

* Run the simulation for a given element number :math:`n`, e.g. :math:`n=1`

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/poisson_3d/build_release
        
        # case (a) - cube shaped elements
        ./linear_regular ../settings_poisson_3d.py 1
        ./quadratic_regular ../settings_poisson_3d.py 1
        
        # case (b) - cuboid shaped elements
        ./linear_structured ../settings_poisson_3d.py 1
        ./quadratic_structured ../settings_poisson_3d.py 1

* Compute the error for the current simulation results.

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/poisson_3d
        python3 compute_error.py

2D Diffusion
----------------------


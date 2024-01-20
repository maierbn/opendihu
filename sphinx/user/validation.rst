Validation
=======================

On this page, we present numeric experiments for problems for which analytic solutions are known. By such comparison, we can validate the implemented solvers in OpenDiHu.

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

Next, the parabolic diffusion equation with a constant coefficient is used to validate the respective solvers.

Problem definition
^^^^^^^^^^^^^^^^^^^^^^

The 2D Diffusion problem is formulated as

.. math::

  -\dfrac{\partial}{\partial t}u(\textbf{x},t) + D\,\Delta u(\textbf{x},t) = 0 \quad \text{on }\textbf{x} \in \Omega = [0, 10] \times [0,10].

The diffusion factor is chosen to be constant :math:`D=3`.

The initial values are given by the following function,

.. math::

    u(\textbf{x},0) = 1 + \textrm{cos}(\textrm{min}(\vert\textbf{x} - \textbf{p}\vert,\pi)), \quad \textbf{p} = (2,5)^\top.

This function describes a peak of height 2 at :math:`\textbf{p}=(2,5)^\top`. The function decays in radial direction away from :math:`\textbf{p}` and reaches constant zero at the boundary :math:`\partial \Omega`.

We consider two scenarios \(a\) and \(b\). In \(a\), we impose homogeneous Dirichlet boundary conditions at :math:`x=0`,

.. math::

    u(x=0,y,t) = 0, \quad (\text{where }\textbf{x} = (x,y)^\top).

In scenario \(b\) we specify homogeneous Neumann boundary conditions at :math:`x=0`, 

.. math::

    \dfrac{\partial}{\partial x} u(x=0,y,t) = 0.

Analytic solution
^^^^^^^^^^^^^^^^^^^^^^

It is known that a solution to the governing diffusion equation on the 2D infinite domain :math:`\Omega_\infty = \mathbb{R}^2` is given by [Ursell2016]_,

.. math::
    :label: diffusion_eq_solution

    u(\textbf{x},t) = \displaystyle\int G(\textbf{x}, \textbf{x}', t)\,u(\textbf{x}', 0)\,\mathrm{d}x,

with *Green's Function*

.. math::

    G(\textbf{x}, \textbf{x}', t) = \dfrac{\mathrm{exp}\Big(-\frac{|\textbf{x}-\textbf{x}'|^2}{4\,D\,t}\Big)}{4\,\pi\,D\,t}.

For scenario \(a\) with homogeneous Dirichlet boundary at :math:`x=0`, we can construct a new function :math:`G_D` as sum of :math:`G(x)` 
and its mirrored counterpart :math:`-G(-x)` whose graph is mirrored around the :math:`x=0` axis,

.. math::

    G_D(x,y,x',y',t) = G(x,y,x',y',t) - G(x,y,-x',y',t)

This function is zero for :math:`x=0`.

Similarly, for the Neumann boundary in scenario \(b\), we construct a function 

.. math::

    G_N(x,y,x',y',t) = G(x,y,x',y',t) + G_D(x,y,-x',y',t).

We get an analytic solution for scenario \(a\) by replacing :math:`G` by :math:`G_D` in Eq. :eq:`diffusion_eq_solution` and for scenario \(b\) by replacing :math:`G` by :math:`G_N` in :eq:`diffusion_eq_solution`.
This solution, however, is correct only for a problem on :math:`\Omega = [0, \infty) \times (-\infty, \infty)`.


Numeric solution 
^^^^^^^^^^^^^^^^^^^^^

With our finite element solver, we can only discretize a finite domain. 
The considered problem on :math:`\Omega = [0,10]\times [0,10]` with the given initial values is specified such that large changes in function value :math:`u` 
mainly occur in the interior of the domain, away from the boundary at :math:`x=10, y=0,` and :math:`y=10`.

:Numref:`diffusiondirichlet` and :numref:`diffusionneumann` show the initial values :math:`u(\textbf{x}, 0)` and the values for :math:`t=1`, 
:math:`u(\textbf{x}, 1)` for scenarios \(a\) and \(b\), respectively.

.. _diffusiondirichlet:
.. figure:: /user/validation/diffusion_dirichlet.png
  :width: 50%

  Scenario \(a\) with homogeneous Dirichlet boundary conditions at :math:`x=0`, for :math:`t=0` (top) and :math:`t=1` (bottom).

.. _diffusionneumann:
.. figure:: /user/validation/diffusion_neumann.png
  :width: 50%

  Scenario \(b\) with homogeneous Neumann boundary conditions at :math:`x=0`, for :math:`t=0` (top) and :math:`t=1` (bottom).

It can also be seen how the different boundary conditions affect the solution value. While the Dirichlet boundary conditions in scenarios \(a\) absorbs the concentration :math:`u` at :math:`x=0`, 
the homogeneous Neumann boundary condition serves as a *reflection boundary* which leads to a constant concentration gradient orthogonal to the wall.
The comparison of the state at :math:`t=0.5` in the following images (Scenario \(a\) left, scenario \(b\) right) shows that concentration decreases faster with the Dirichlet boundary condition.

.. |diffusion_dirichlet_0_5| image:: /user/validation/diffusion_dirichlet_0_5.png
    :width: 48%

.. |diffusion_neumann_0_5| image:: /user/validation/diffusion_neumann_0_5.png
    :width: 48%

|diffusion_dirichlet_0_5| |diffusion_neumann_0_5|

We disretize the problem in space by :math:`10\times 10` quadratic Finite Elements and in time using the Crank-Nicolson scheme with time step width :math:`dt = 1e-5`. 
The linear system of equations in every timestep is solved by a geometric multi-grid solver.

Comparison of numeric and analytic solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We compare the solution value :math:`u(\textbf{x},t)` over time at selected nodes of the Finite Element mesh. 
The quadratic mesh has :math:`21 \times 21` nodes that span the :math:`[0,10] \times [0,10]` domain.
:numref:`diffusion_error_dirichlet` and :numref:`diffusion_error_neumann` show the results for scenario \(a\) and scenario \(b\), the selected nodes are given in the legend.

.. _diffusion_error_dirichlet:
.. figure:: /user/validation/diffusion_error_dirichlet.png
  :width: 80%

  Scenario (a): Numeric (solid lines) and analytic solution (dashed lines) for selected nodes in the mesh.

.. _diffusion_error_neumann:
.. figure:: /user/validation/diffusion_error_neumann.png
  :width: 80%

  Scenario (b): Numeric (solid lines) and analytic solution (dashed lines) for selected nodes in the mesh.

It can be seen that the analytic and numeric plots essentially coincide. 
At some points, the solution starts to differ during the end of the simulation time span, which can be explained by the analytic solution being correct for the infinite domain.

This result shows that the diffusion problem solver of OpenDiHu is validated.

How to reproduce
^^^^^^^^^^^^^^^^^^^^^

* Build the example

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/diffusion_2d
        mkorn && sr       # build

* Run the simulations

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/diffusion_2d/build_release
        
        # scenario (a)
        ./diffusion2d_quadratic ../settings_diffusion_2d_dirichlet.py
        
        # scenario (b)
        ./diffusion2d_quadratic ../settings_diffusion_2d_neumann.py

* Visualize the numeric solutions (make sure that the ``$OPENDIHU_HOME/scripts`` directory is in the ``PYTHONPATH`` environment variable)

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/diffusion_2d/build_release
        
        # scenario (a)
        plot out_dirichlet/*

        # scenario (b)
        plot out_neumann/*


* Calculate the analytic solution and create the plot

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/diffusion_2d
        
        # scenario (a)
        python3 compute_error_dirichlet.py

        # scenario (b)
        python3 compute_error_neumann.py

.. [Ursell2016] `Ursell et al. 2016, The Diffusion Equation/A Multi-dimensional Tutorial <https://www.yumpu.com/en/document/read/7921375/the-diffusion-equation-a-multi-dimensional-tutorial-california->`_

CellML models
-----------------------------------

We solve CellML models with `OpenCOR <https://opencor.ws/>`_ and OpenDiHu and compare the results. 
`CellML <https://www.cellml.org/>`_ is a description standard for differential-algebraic systems of equations used in the bioengineering community. The XML based CellML file format can be parsed and solved by both OpenCOR and OpenDiHu.

Hodgkin-Huxley (1952)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most basic electrophysiology model was formulated by `Hodgkin and Huxley (1952) <https://royalsocietypublishing.org/doi/abs/10.1098/rspb.1952.0054>`_. It contains four state variables :math:`V_m`, :math:`n`, :math:`h`, :math:`m`.
We set a constant stimulation current of :math:`i_\text{Stim} = 10 μA/cm^2` and solve the system for :math:`t_\text{end}=35 ms` with Heun's method and a time step witdh of :math:`dt=1e-5 ms`.

As a reference, the same simulation is performed with OpenCOR. The resulting values show good agreement.

.. _hodgkin_huxley_opencor:
.. figure:: /user/validation/hodgkin_huxley_opencor.png
  :width: 70%

  Visualization of the results in the OpenCOR GUI
  
.. _hodgkin_huxley_opendihu:
.. figure:: /user/validation/hodgkin_huxley_opendihu.png
  :width: 70%

  Visualization of the results from OpenDiHu with the OpenDiHu plot script
  
.. _hodgkin_huxley_compared:
.. figure:: /user/validation/hodgkin_huxley_compared.png
  :width: 70%

  Comparison of results computed by OpenCOR (dashed lines) and OpenDiHu (solid lines), which show good agreement.
  
The :math:`\mathcal{L}_2`-errors for the solved variables over the entire timespan are given below.

.. list-table::
  :widths: 5 10
  :header-rows: 1

  * - variable
    - :math:`\mathcal{L}_2`-error
  * - :math:`V_m`
    - 2.0e-03
  * - :math:`n`
    - 3.3e-06
  * - :math:`h`
    - 4.2e-06
  * - :math:`m`
    - 1.9e-05

Shorten, Ocallaghan, Davidson, Soboleva (2007)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another relevant CellML model is the model by `Shorten, Ocallaghan, Davidson, Soboleva (2007) <https://link.springer.com/article/10.1007/s10974-007-9125-6>`_.
It consists of 58 ordinary differential equations and 77 algebraic equations that have to be solved in time. 
We repeat the previous study with this model using the same numerical parameters and present the result in the following.
The stimulation current ``wal_environment/I_HH`` is set to 50.

.. _cellml_shorten_ocallaghan_comparison:
.. figure:: /user/validation/cellml_shorten_ocallaghan_comparison.png
  :width: 100%

  Comparison of results computed by OpenCOR (dashed lines) and OpenDiHu (solid lines), which show good agreement.

The :math:`\mathcal{L}_2`-errors over all timesteps for all state variables are given below. 
It can be seen that they show good agreement with the reference simulation.


.. list-table::
  :header-rows: 1

  * - variable
    - :math:`\mathcal{L}_2`-error
    - variable
    - :math:`\mathcal{L}_2`-error
    - variable
    - :math:`\mathcal{L}_2`-error
    - variable
    - :math:`\mathcal{L}_2`-error
  * - A_1
    - 5.7e-06
    - A_2
    - 3.5e-06
    - ATP1
    - 1.3e-03
    - ATP2
    - 3.6e-04
  * - Ca_1
    - 9.9e-04
    - Ca_2
    - 1.1e-04
    - Ca_ATP1
    - 1.4e-03
    - Ca_ATP2
    - 3.7e-04
  * - Ca_CaT2
    - 2.5e+01
    - Ca_Cs1
    - 2.6e-02
    - Ca_Cs2
    - 2.7e-02
    - Ca_P1
    - 2.9e-03
  * - Ca_P2
    - 3.3e-04
    - Ca_SR1
    - 1.8e-02
    - Ca_SR2
    - 1.8e-02
    - Ca_T_2
    - 2.5e+01
  * - D_0
    - 9.2e-06
    - D_1
    - 1.9e-05
    - D_2
    - 5.7e-05
    - Mg1
    - 2.9e-03
  * - Mg2
    - 2.8e-03
    - Mg_ATP1
    - 2.9e-03
    - Mg_ATP2
    - 2.8e-03
    - Mg_P1
    - 2.8e-04
  * - Mg_P2
    - 2.9e-04
    - P
    - 2.9e-06
    - P_C_SR
    - 2.9e-07
    - P_SR
    - 2.9e-08
  * - x_1
    - 0.0e+00
    - x_2
    - 0.0e+00
    - C_0
    - 1.9e-05
    - C_1
    - 1.3e-05
  * - C_2
    - 5.9e-06
    - C_3
    - 8.1e-07
    - C_4
    - 3.2e-08
    - O_0
    - 3.8e-11
  * - O_1
    - 6.3e-10
    - O_2
    - 7.3e-09
    - O_3
    - 2.5e-08
    - O_4
    - 2.2e-08
  * - K_e
    - 3.0e-06
    - K_i
    - 2.8e-04
    - K_t
    - 1.3e-05
    - Na_e
    - 3.2e-04
  * - Na_i
    - 2.7e-05
    - Na_t
    - 2.8e-04
    - h_K
    - 3.3e-07
    - n
    - 1.6e-05
  * - h
    - 8.2e-06
    - m
    - 8.1e-05
    - S
    - 3.0e-07
    - h_K_t
    - 6.1e-05
  * - n_t
    - 1.5e-05
    - h_t
    - 1.3e-05
    - m_t
    - 8.2e-05
    - S_t
    - 2.8e-07
  * - vS
    - 5.2e-03
    - vT
    - 4.4e-03
    -
    -
    -
    -

How to reproduce
^^^^^^^^^^^^^^^^^^^^^

* Build the example

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/cellml
        mkorn && sr       # build

* Run the simulations

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/cellml/build_release
        
        ./hodgkin_huxley ../settings_cellml.py
        ./shorten ../settings_cellml.py

* Visualize the numeric solutions (make sure that the ``$OPENDIHU_HOME/scripts`` directory is in the ``PYTHONPATH`` environment variable)

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/cellml/build_release
        
        plot out_hodgkin_huxley/*
        plot out_shorten/*


* Compute the L2 errors and generate the comparison plot

    The reference values were simulated in OpenCOR (set ``wal_environment/I_HH`` manually to 50) and then exported as csv files, only selected the state variables to be included in the file.
    The repository contains these reference values in the directory ``$OPENDIHU_HOME/examples/validation/cellml/`` in the files ``hodgkin_huxley_1952_data_istim10_opencor.csv`` and ``shorten_ocallaghan_davidson_soboleva_2007_no_stim_data_ihh50_opencor.csv``.

    .. code-block:: bash

        cd $OPENDIHU_HOME/examples/validation/cellml
        
        python3 compute_error_hodgkin_huxley.py
        python3 compute_error_shorten.py

Nonlinear solid mechanics 
-----------------------------

To validate the implementation of the nonlinear incompressible solid mechanics solver, a comparison against results from the `FEBio <https://febio.org/>`_ solver are conducted.
Two test case of a tensile test and a shear test are simulated with both solvers.
Three different formulations of incompressible hyperelasticity in OpenDiHu are compared with the reference solution from FEBio and found to be equal.

The description of the setup and the results can be found in the `dissertation <https://arxiv.org/abs/2107.07104>`_, chapter 8.2.2 [Maier2021]_.

.. _solid_mechanics_validation_tensile_test:
.. figure:: /user/validation/solid_mechanics_validation_tensile_test.png
  :width: 100%

  Resulting stresses in the tensile test showing good match between the different formulations with OpenDiHu and the reference solver FEBio
  

.. _solid_mechanics_shear_test:
.. figure:: /user/validation/solid_mechanics_shear_test.png
  :width: 70%

  Resulting stress components in the shear test with good match to the reference solver FEBio
  


.. [Maier2021] `Maier (2021), Scalable Biophysical Simulations of the Neuromuscular System <https://arxiv.org/abs/2107.07104>`_, Ph.D. thesis, University of Stuttgart, 2021

How to reproduce
^^^^^^^^^^^^^^^^^^^^^

* Build the example

    .. code-block:: bash

        # tensile test
        cd $OPENDIHU_HOME/opendihu/examples/solid_mechanics/tensile_test
        mkorn && sr

        # shear test
        cd $OPENDIHU_HOME/opendihu/examples/solid_mechanics/shear_test
        mkorn && sr

* Run the simulations

    For the FEBio simulations to work you need FEBio 3 installed (such that it can be called from the command line as ``febio3``). If FEBio is not installed, the OpenDiHu simulations will still run, but the reference data will not be available in the plot.

    .. code-block:: bash

        # tensile test
        cd $OPENDIHU_HOME/examples/solid_mechanics/tensile_test/build_release
        ../run_force.sh

        # shear test
        cd $OPENDIHU_HOME/examples/solid_mechanics/shear_test/build_release
        ../run_force.sh

* Visualize the results

    .. code-block:: bash

        # tensile test
        cd $OPENDIHU_HOME/examples/solid_mechanics/tensile_test
        python3 plot_validation.py

        # shear test
        cd $OPENDIHU_HOME/examples/solid_mechanics/shear_test
        python3 plot_validation.py


Excited cuboid muscle
-----------------------

For more complex multi-scale models, no analytic solutions are available and a direct verification of the solver becomes impossible.
A qualitative validation can be carried out by comparing to other groups simulation. In the following, we compute a scenario similar to the one in [Heidlauf2016]_, chapter 7.3.1.

A cuboid muscle of size :math:`2.9 \times 1.2 \times 6` cm has :math:`31 \times 13` embedded muscle fibers.
The fibers are organized in 10 motor units which are activated from a motor neuron model, as described in [Heidlauf2016]_. We use the same discharge times as the author.
We the discretize the fibers by 145 linear Finite Elements as in [Heidlauf2016]_. No information is given about the 3D domain discretization, for which we use :math:`15 \times 5 \times 33 = 2475` quadratic Finite Elements.
A difference between our simulation and the one in [Heidlauf2016]_ is that we omit the fat layer and use a more realistic dynamic solid mechanics formulation instead of the quasi-static formulation.

The contraction state in our simulation is shown below and can be compared to Fig. 7.4 in [Heidlauf2016]_. (Note that we set the `"displacementsScalingFactor"` to 5, equal to the visualization in the referenced thesis.)
The contraction shows qualitative agreement.

.. [Heidlauf2016] `Heidlauf (2016) Chemo-electro-mechanical modelling of the neuromuscular system <http://dx.doi.org/10.18419/opus-658>`_, Ph.D. thesis, University of Stuttgart

.. |fibers004| image:: /user/validation/fibers004.png
    :width: 70%

.. |fibers020| image:: /user/validation/fibers020.png
    :width: 70%

.. |fibers040| image:: /user/validation/fibers040.png
    :width: 70%

.. |fibers080| image:: /user/validation/fibers080.png
    :width: 70%

.. |fibers120| image:: /user/validation/fibers120.png
    :width: 70%

.. |fibers200| image:: /user/validation/fibers200.png
    :width: 70%

|fibers004| |fibers020|
|fibers040| |fibers080|
|fibers120| |fibers200|

:numref:`fibers_stretch` evaluates the amount of contraction that the muscle undergoes in this scenario. The orange lines shows the average over all nodes of the 3D mesh, the light orange range corresponds to the 25\%-75\% quartile range of values.
At 240ms, the average stretch is circa 90% which means that the muscle contracted by 10% in this time.

.. _fibers_stretch:
.. figure:: /user/validation/fibers_stretch.png
  :width: 100%

  Stretch over time of the muscle tissue. 


More details about the scenario in OpenDiHu can be found here:

.. code-block::

    0/16 : This is opendihu 1.3, built Aug 29 2023, C++ 201402, GCC 9.4.0, current time: 2023/8/30 11:10:01, hostname: pcsgs05, n ranks: 16      
    0/16 : Open MPI v4.0.3, package: Debian OpenMPI, ident: 4.0.3, repo rev: v4.0.3, Mar 03, 2020                                                
    0/16 : File "../settings_cuboid_muscle.py" loaded.                                                                                           
    0/16 : ---------------------------------------- begin python output ----------------------------------------                                 
    Loading variables from "heidlauf.py".                                                                                                        
    scenario_name: heidlauf,  n_subdomains: 2 2 4,  n_ranks: 16,  end_time: 4000.0                                                               
    dt_0D:           1e-04, diffusion_solver_type:      cg                                                                                       
    dt_1D:           1e-04, potential_flow_solver_type: gmres                                                                                    
    dt_splitting:    1e-04, emg_solver_type:            cg, emg_initial_guess_nonzero: False                                                     
    dt_3D:           1e+00, paraview_output: True                                                                                                
    output_timestep: 1e+00  stimulation_frequency: 1.0 1/ms = 1000.0 Hz                                                                          
    fast_monodomain_solver_optimizations: True, use_analytic_jacobian: True, use_vc: True                                                        
    fiber_file:              cuboid.bin                                                                                                          
    fat_mesh_file:           ../../../../input/13x13fibers.bin_fat.bin                                                                           
    cellml_file:             ../../../electrophysiology/input/new_slow_TK_2014_12_08.c                                                           
    fiber_distribution_file: ../../../electrophysiology/input/MU_fibre_distribution_10MUs.txt                                                    
    firing_times_file:       ../../../electrophysiology/input/MU_firing_times_heidlauf_10MU.txt                                                  
    ********************************************************************************                                                             
    prefactor: sigma_eff/(Am*Cm) = 0.0132 = 3.828 / (500.0*0.58)                                                                                 
    create cuboid.bin with size [2.9,1.2,6], n points [31,13,145]                                                                                
    diffusion solver type: cg                                                                                                                    
    n fibers:              403 (31 x 13), sampled by stride 2 x 2                                                                                
    n points per fiber:    145, sampled by stride 4                                                                                              
    16 ranks, partitioning: x2 x y2 x z4                                                                                                         
    31 x 13 = 403 fibers, per partition: 14 x 6 = 84                                                                                             
    per fiber: 1D mesh    nodes global: 145, local: 36                                                                                           
      sampling 3D mesh with stride 2 x 2 x 4                                                                                                     
      distribute_nodes_equally: True                                                                                                             
    quadratic 3D mesh    nodes global: 15 x 5 x 33 = 2475, local: 8 x 2 x 8 = 128                                                               
    quadratic 3D mesh elements global: 7 x 2 x 16 = 224, local: 4 x 1 x 4 = 16                                                                  
    number of degrees of freedom:                                                                                                                
                        1D fiber:        145  (per process: 36)                                                                                  
                0D-1D monodomain:       8120  (per process: 2016)                                                                                
    all fibers 0D-1D monodomain:    3272360  (per process: 169344)                                                                              
                    3D bidomain:       2475  (per process: 128)                                                                                 
                          total:    3274835  (per process: 169472)                                                                              
    Debugging output about fiber firing: Taking input from file "../../../electrophysiology/input/MU_firing_times_heidlauf_10MU.txt"             
    First stimulation times                                                                                                                      
        Time  MU fibers                                                                                                                          
        0.00   2 [60, 80, 81, 101, 104, 118, 127, 191, 228, 275] (only showing first 10, 19 total)                                               
        2.00   3 [1, 12, 17, 40, 75, 76, 82, 94, 98, 108] (only showing first 10, 22 total)                                                      
        3.00   5 [7, 10, 19, 21, 24, 32, 35, 43, 62, 68] (only showing first 10, 42 total)                                                       
        7.01   4 [14, 16, 27, 39, 52, 106, 111, 116, 125, 146] (only showing first 10, 35 total)                                                 
        8.01   7 [0, 2, 9, 15, 18, 28, 33, 36, 45, 55] (only showing first 10, 53 total)                                                         
        9.01   1 [29, 49, 58, 102, 131, 134, 153, 182, 184, 232] (only showing first 10, 14 total)                                               
      16.02   6 [26, 31, 44, 51, 53, 56, 57, 65, 67, 83] (only showing first 10, 46 total)                                                      
      35.04   8 [3, 13, 23, 25, 30, 38, 41, 54, 71, 78] (only showing first 10, 69 total)                                                       
      never stimulated: MU   0, fibers [37, 46, 47, 48, 69, 74, 89, 93, 138, 213] (only showing first 10, 23 total)                              
    stimulated MUs: 8, not stimulated MUs: 1                                                                                 

How to reproduce
^^^^^^^^^^^^^^^^^^^^^

* Build the example

    .. code-block:: bash

        cd $OPENDIHU_HOME/opendihu/examples/validation/cuboid_muscle/
        mkorn && sr

* Run the simulation

    .. code-block:: bash

        cd $OPENDIHU_HOME/opendihu/examples/validation/cuboid_muscle/build_release
        mpirun -n 16 ./cuboid_muscle ../settings_cuboid_muscle.py heidlauf.py


.. Cuboid muscle with EMG
.. -------------------------


.. How to reproduce
.. ^^^^^^^^^^^^^^^^^^^^^

.. * Build the example

..     .. code-block:: bash
        
..         cd $OPENDIHU_HOME/opendihu/examples/electrophysiology/fibers/fibers_fat_emg_contraction/
..         mkorn && sr

.. * Run the simulation

..     .. code-block:: bash

..         cd $OPENDIHU_HOME/opendihu/examples/electrophysiology/fibers/fibers_fat_emg_contraction/build_release
..         ./fibers_fat_emg_contraction ../settings_fibers_fat_emg_contraction.py fibers.py

.. .. code-block::

..     Open MPI v4.0.3, package: Debian OpenMPI, ident: 4.0.3, repo rev: v4.0.3, Mar 03, 2020
..     File "../settings_fibers_fat_emg_contraction.py" loaded.
..     ---------------------------------------- begin python output ----------------------------------------
..     Loading variables from "fibers.py".
..     scenario_name: fibers,  n_subdomains: 1 1 1,  n_ranks: 1,  end_time: 5000.0
..     dt_0D:           2.5e-05, diffusion_solver_type:      cg
..     dt_1D:           2.5e-05, potential_flow_solver_type: gmres, approx. exp.: True
..     dt_splitting:    2.5e-05, emg_solver_type:            cg, emg_initial_guess_nonzero: False
..     dt_3D:           1.0e-03, paraview_output: True, optimization_type: vc
..     dt_elasticity:   1e-01    elasticity solver: lu, preconditioner: none
..     fiber_file:              cuboid.bin
..     fat_mesh_file:           cuboid_fat2.bin
..     cellml_file:             /data/scratch/maierbn/opendihu/examples/electrophysiology/input/new_slow_TK_2014_12_08.c
..     fiber_distribution_file: /data/scratch/maierbn/opendihu/examples/electrophysiology/input/MU_fibre_distribution_10MUs.txt
..     firing_times_file:       /data/scratch/maierbn/opendihu/examples/electrophysiology/input/MU_firing_times_heidlauf_10MU.txt
..     ********************************************************************************
..     n fibers:              403 (31 x 13), sampled by stride 1 x 1
..     n points per fiber:    145, sampled by stride 10
..     1 rank, partitioning: x1 x y1 x z1
..     31 x 13 = 403 fibers, per partition: 30 x 12 = 360
..     per fiber: 1D mesh    nodes global: 145, local: 145
..       sampling 3D mesh with stride 1 x 1 x 10 
..         linear 3D mesh    nodes global: 31 x 13 x 15 = 6045, local: 31 x 13 x 15 = 6045
..         linear 3D mesh elements global: 30 x 12 x 14 = 5040, local: 30 x 12 x 14 = 5040
..     quadratic 3D mesh    nodes global: 31 x 13 x 15 = 6045, local: 31 x 13 x 15 = 6045
..     quadratic 3D mesh elements global: 15 x 6 x 7 = 630, local: 15 x 6 x 7 = 630
..     number of degrees of freedom:
..                         1D fiber:        145  (per process: 145)
..                 0D-1D monodomain:       8120  (per process: 8120)
..     all fibers 0D-1D monodomain:    3272360  (per process: 2923200)
..                     3D bidomain:       6045  (per process: 6045)
..                           total:    3278405  (per process: 2929245)
..         fat mesh, n points total:    1935 (43 x 3 x 15), (per process: 43 x 3 x 15 = 1935)
..       sub-sampling 3D elasticity mesh with factors 0.7, 0.7, 0.7 
..       elasticity quadratic 3D meshes:
..       muscle:             nodes global: 21 x 9 x 11 = 2079, local: 21 x 9 x 11 = 2079
..               quadratic elements global: 10 x 4 x 5 = 200, local: 14 x 1 x 5 = 70
..       fat and skin layer: nodes global: 29 x 3 x 11 = 957, local: 29 x 3 x 11 = 957
..               quadratic elements global: 14 x 1 x 5 = 70, local: 14 x 1 x 5 = 70
..     Python config parsed in 0.2s.
..     ----------------------------------------- end python output -----------------------------------------
..     Read from file "cuboid.bin", 16651 collective chunks.
..     done.
..     Note, setting regularization to [tol,eps] = [1e-2,1e-1]. Add `"regularization": None` in the settings to disable.
..     Initialize 1 global instances (1 local).
..     CellML file "/data/scratch/maierbn/opendihu/examples/electrophysiology/input/new_slow_TK_2014_12_08.c" with 57 states, 71 algebraics, specified 2 parameters: 
..       parameter 0 maps to "wal_environment/I_HH" (CONSTANTS[54]), initial value: 0, 
..       parameter 1 maps to "razumova/L_S" (CONSTANTS[67]), initial value: 1

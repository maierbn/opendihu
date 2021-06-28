
Implement New Physics
========================

This page is a tutorial about how to implement new a physics model that is not yet possible with the C++ core.
It covers the setup steps to create a new solver class, how to use the existing integration functionality, how to use the available infrastructure such as output writers, Petsc vectors etc. and how to create an example that uses the new model.
This will be done by considering an exemplary *advection-diffusion equation*.

A "reference solution" was created, see the ``diffusion_advection_solver.tpp`` files under ``specialized_solvers/darcy`` and the example under ``examples/darcy/diffusion_advection2d``.

1. *Complete the formulation.* 

  The first step is to derive the formulation and be clear about what exactly the algorithm has to do in order to solve the model. In this example case, we want to solve the following equation.
    
  .. math::
    \frac{∂u(x,t)}{∂t} = Δu(x,t) - ∇\cdot (vu) + f(x,t),
    
    
  where :math:`u \in \mathbb{R}` is the scalar solution variable, :math:`v \in \mathbb{R}^d` is a prescribed velocity field and :math:`f(x,t) \in \mathbb{R}` is a reaction term.
  The dimensionality, :math:`d` of the domain :math:`\Omega \subset \mathbb{R}^d \ni x`, should be any of 1,2 or 3.

  We want to use the Finite Element method, therefore we derive the weak form. For simplicity, we discretize in time using the forward Euler scheme.

  .. math::
    \int_\Omega u_t \, \phi \, dx &= -\int_\Omega ∇u \cdot ∇\phi \,dx + \int_\Omega (vu) \cdot ∇\phi \, dx + \int_\Omega f \, \phi \, dx, \quad \forall \phi \in H^{1}_0(\Omega) \\
    \Leftrightarrow \quad &\frac1{dt} \sum_{i=1}^N (u_i^{(t+1)} - u_i^{(t)})\,\int_\Omega \phi_i \, \phi_j \,dx \\
    \quad &= - \sum_{i=1}^N u_i^{(t)} \int_\Omega (∇\phi_i \cdot ∇\phi_j + (v\phi_i)\cdot ∇\phi_j)\, dx + \int_\Omega f \,\phi_j \,dx \quad \forall j

  By using matrix notation with :math:`K,M` and :math:`V` for the stiffness, mass and velocity matrices and the vector of degrees of freedom, :math:`u`, we get

  .. math::
    \frac1{dt} M (u^{(t+1)} - u^{(t)}) &= K\,u^{(t)} - V\,u^{(t)} - M\,f\\
    \Leftrightarrow \quad u^{(t+1)} &= u^{(t)} + dt\cdot( M^{-1}Ku^{(t)} - M^{-1}Vu^{(t)} - f),
    
  where

  .. math::
    M_{ij} = \int_\Omega \phi_i\,\phi_j \,dx, \quad K_{ij} = -\int_\Omega ∇\phi_i \cdot ∇\phi_j \,dx, \quad
    V_{ij} = -\int_\Omega \phi_i\,(v \cdot ∇\phi_j) \,dx.
  
2. *Create the solver class.* 

  Next, we decide how to name the new solver class. Here, we choose `AdvectionDiffusionSolver`. There are two classes that can be used as template to start the new solver from. 

  There is ``specialized_solver/my_new_static_solver.h`` for static problem and ``specialized_solver/my_new_timestepping_solver.h`` for transient problem, i.e. timestepping schemes.
  In our case we need the first, so we copy ``specialized_solver/my_new_timestepping_solver.h`` to ``specialized_solver/advection_diffusion_solver.h`` and ``specialized_solver/my_new_timestepping_solver.tpp`` to ``specialized_solver/advection_diffusion_solver.tpp``.
  
  Replace all occurences of ``MyNewTimesteppingSolver`` with ``AdvectionDiffusionSolver``. Read the comments to understand what this template already contains and what still needs to be done.

  We want the reaction term, :math:`f` to be specified by a python callback function. This functionality is implemented by the `PrescribedValues` class which we want to reuse. Such a class should be a nested `solver` of the `AdvectionDiffusionSolver`. 
  Therefore, we need to implement one template type. Replace all ``template<class TimeStepping>`` by ``template<class PrescribedValues>`` and all ``AdvectionDiffusionSolver<TimeStepping>::`` by ``AdvectionDiffusionSolver<PrescribedValues>::``. 
  Note this is only a renaming but it would be confusing to call the nested `PrescribedValues` instead `class TimeStepping`.

3. *Prepare the data class.* 

  Add scalar and vector-valued field variables for the reaction term, :math:`f \in \mathbb{R}^1` and the velocity field, :math:`v \in \mathbb{R}^d`.

4. *Implement the assembly of matrices.*

  There is already code to do the quadrature of the integral terms. Open the file `spatial_discretization/finite_element_method/01_stiffness_matrix_integrate.tpp` and copy the `setStiffnessMatrix()` method to the `AdvectionDiffusionSolver`.
  Rename it, e.g. ``setSystemMatrix``. This method computes only the stiffness matrix, :math:`K`. Carefully read the code and try to understand it. Change it to produce the velocity matrix, :math:`V`, instead.

  In order to get the other matrices, :math:`K` and :math:`M^{-1}`, you can use an object of type ``FiniteElementMethod`` as member variable of ``AdvectionDiffusionSolver``. This class already implements assembly of the normal stiffness matrix, as well as an inverse lumped mass matrix, to be used for :math:`M^{-1}`.
  Inspect the data object of the ``FiniteElementMethod`` which is in `data_management/finite_element_method/finite_elements_base.h`. You can see that there are getter methods `stiffnessMatrix()` and `inverseLumpedMassMatrix()`. In order for them to be available, the initialization method `initializeInverseLumpedMassMatrix` has to be called.

5. *Create an example and test it.*

6. *Create a unit test.*  

  The last step is to create a unit test to preserve the functionality that is working now also when the code is further changed. The unit test can be a regression test that compares the output to a reference output that is now created.

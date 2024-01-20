Connector Slots
===================================

Each solver has multiple *connector slots* where field variables are exposed to surrounding solvers. 
When using a :doc:`/settings/coupling` or :doc:`/settings/splitting` scheme, data will be transferred between the two involved solvers.

After solution of term 1 is complete, all slots of term 2 will get the values of the connected slots from term 1. After the solution of term 2 is complete, the opposite connections will be used.

The connections between the slots of two solvers have to be specified in the settings. There are two possibilities

* Either in the Coupling or Splitting scheme itself using the ``connectedSlotsTerm1To2`` and ``connectedSlotsTerm2To1`` options.
* Or globally with the ``"connectedSlots"`` option and referencing the slots by slot names.

The first possibility is shown below.
The following is an example of two solvers within a `StrangSplitting`:

.. code-block:: python

  "solverStructureDiagramFile": "solver_structure.txt",     # filename of file that will contain a visualization of the solver structure and data mapping  
  "StrangSplitting": {
    #"timeStepWidth" and other options
    
    "connectedSlotsTerm1To2": [1,0],         # list of slots of term 2 that are connected to the slots of term 1
    "connectedSlotsTerm2To1": [0,None,1],    # list of slots of term 1 that are connected to the slots of term 2, None means slot is not connected
    
    #"connectedSlotsTerm1To2": {0:1, 1:0},   # alternative form
    #"connectedSlotsTerm2To1": {0:0, 2:1},
    
    "Term1": { 
      # settings of first solver
    },
    "Term2": { 
      # settings of second solver
    }
  }

`connectedSlotsTerm1To2` and `connectedSlotsTerm2To1`
--------------------------------------------------------
The options ``connectedSlotsTerm1To2`` and ``connectedSlotsTerm2To1`` specify which slots should be connected to each other. The first specifies to which slots of term 2 the slots of term 1 will be connected to. The second vice-versa.

They can either be given as a list or as a dict. If a list is given, the items in the list are the slots of the other solver and correspond to the own slots 0,1,2... Slots that should not be connected have to be given as ``None``. But normally, all slots should be connected.

If a dict is given, the mappings can be specified directly, i.e. ``2:1`` means connect the own slot no. 2 to the slot 1 of the other solver.

If two slots of both solvers are connected both ways, the field variables can be reused. Otherwise, the data has to be copied, which is slower. Therefore, it is often beneficial to connect slots both ways even if the reverse direction wouldn't be needed for the data flow, but then copy can be avoided.

There are some special rules, when data will still be copied even if the slots are connected both ways, e.g. if there is a :doc:`/settings/cellml_adapter` and two different states are mapped. This has to do with the internal data representation using the PETSc Vec's.

If copy is avoided can be seen from the *solverStructureDiagramFile* which will be explained below.


solverStructureDiagramFile
--------------------------------------------------------
(string) Filename of a file that will contain a visualization of the output slots connections. This file will be written when the program completes. If the program crashes, the file will also be written, however, it may not contain the full information if not all solvers have been visited prior to the crash.

The following is an example for such a file. It is from the ``fibers_emg`` example. (It is an older version, the recent version looks even fancier) It can be seen that the ``solution.membrane/V`` variable is shared between the CellMLAdapter and the diffusion solver (no copy). Outside this variable is also copied to the StaticBidomainSolver.

.. code-block:: bash

  Solver structure: 

  ├── Coupling                                     
  │  output slots:                                 
  │  solution.membrane/V                     ── ¤0 x
  │                                                
  │  slot connections:                             
  │  0¤ -> ¤0                                      
  │                                                
  │ ├── FastMonodomainSolver                       
  │ │   ("Term1")                                  
  │ │ └── MultipleInstances                        
  │ │ │ ├── StrangSplitting                        
  │ │ │ │  output slots:                           
  │ │ │ │  solution.membrane/V               ── ¤0────┐
  │ │ │ │                                             │
  │ │ │ │  slot connections:                          │
  │ │ │ │  0¤ <=> ¤0                                  │
  │ │ │ │                                             │
  │ │ │ │ ├── MultipleInstances                       │
  │ │ │ │ │   ("Term1")                               │
  │ │ │ │ │ ├── Heun                                  │
  │ │ │ │ │ │  output slots:                          │
  │ │ │ │ │ │  solution.membrane/V           ── ¤0══╗ │
  │ │ │ │ │ │                                       ║ │
  │ │ │ │ │ │ └── CellmlAdapter                     ║ │
  │ │ │ │ │ └                                       ║ │
  │ │ │ │ └                                         ║ │
  │ │ │ │                                           ║ │
  │ │ │ │ └── MultipleInstances                     ║ │
  │ │ │ │ │   ("Term2")                             ║ │
  │ │ │ │ │ ├── ImplicitEuler                       ║ │
  │ │ │ │ │ │  output slots:                        ║ │
  │ │ │ │ │ │  solution.0                    ── ¤0══╝ │
  │ │ │ │ │ │                                         │
  │ │ │ │ │ │ ├── FiniteElementMethod                 │
  │ │ │ │ │ │ │  output slots:                        │
  │ │ │ │ │ │ │  solution.0                  ── ¤0 x  │
  │ │ │ │ │ │ │                                       │
  │ │ │ │ │ └                                         │
  │ │ │ │ └                                           │
  │ │ │ └                                             │
  │ │ └                                               │
  │ └                                                 │
  │                                                   │
  │ ├── StaticBidomainSolver                          │
  │ │   ("Term2")                                     │
  │ │  output slots:                                  │
  │ │  Vm.0                                  ── ¤0<───┘
  │ │                                              
  │ │ ├── FiniteElementMethod                      
  │ │ │   ("PotentialFlow")                        
  │ │ │  output slots:                             
  │ │ │  solution.0                          ── ¤0 x
  │ │ │                                            
  │ │                                              
  │ │ ├── FiniteElementMethod                      
  │ │ │   ("Activation Transmembrane")             
  │ │ │  output slots:                             
  │ │ │  solution.0                          ── ¤0 x
  │ │ │                                            
  │ │                                              
  │ │ ├── FiniteElementMethod                      
  │ │ │   ("Activation Extracellular")             
  │ │ │  output slots:                             
  │ │ │  solution.0                          ── ¤0 x
  │ │ │                                            
  │ └                                              
  └                                                
                                                   
  connection types:
    ═══ ... reuse field variable, no copy
    ──> ... copy data in direction of arrow

Using global slot names
-----------------------------------
Another possibility that is advantageous for more complex examples is to specify all slot connections globally.
This required that all connector slots have names assigned. These names have to be set by options in the solvers, usually ``slotNames`` (e.g. :doc:`static_bidomain_solver`, :doc:`muscle_contraction_solver`, :doc:`quasi_static_linear_elasticity_solver`) or ``additionalSlotNames`` (e.g. any :doc:`timestepping_schemes_ode`, :doc:`map_dofs`). For the :doc:`cellml_adapter`, the slot names are directly given in the ``mappings`` option.

Then you can define the option

.. code-block:: python

  config = {
    ...
    "connectedSlots": [
      ("mn_out", "mn"),
      ("in_g",   "in_in"),
      ("msin_s", "msin_i"),
      ("msin_i", "msin_m"),
      ("gt",     "gt_in"),
      ("ms",     "ms_in"),
    ],
    ...
  }

It is a list of tuples with ``("fromName", "toName")`` entries. 

Note that the slot names must be 10 characters long or less. This restriction is because of the solver structure visualization. (Actually they can be any length but only the first 10 characters will be shown in the solver structure file.)

Connection by choosing equal slot names
---------------------------------------------

Slots that have the same name are automatically connected. Thus, it is often enough to just name all slots properly.


Example
------------

The ``examples/electrophysiology/monodomain/hodgkin_huxley`` example yields the following solver structure (file ``solver_structure.txt``):

.. code-block:: bash

  The following data slot connection were given by the setting "connectedSlots":
         h ¤ <─> ¤ h_gate

  The following data slots were connected because the names appeared in both terms of a coupling or splitting scheme:
    m_gate ¤ <─> ¤ m_gate

  Solver structure: 

  ├── StrangSplitting                                                
  │  data slots:                                                     
  │  [a] solution.membrane/V                     ├─────────────── ¤0 x
  │  [a] solution.sodium_channel_m_gate/m        :├────────m_gate ¤1 x
  │  [a] solution.sodium_channel_h_gate/h        ::├───────h_gate ¤2 x
  │  [a] solution.potassium_channel_n_gate/n     :::├──────────── ¤3 x
  │  [a] additionalFieldVariable0                ::::├──────── aa ¤4 x
  │  [a] additionalFieldVariable1                :::::├─────── bb ¤5 x
  │  [a] leakage_current/i_L                     ::::::├───────── ¤6 x
  │  [a] solution                                :::::::├───── vm ¤7 x
  │  [a] additionalFieldVariable0                ::::::::├─m_gate ¤8 x
  │  [a] additionalFieldVariable1                :::::::::├──── h ¤9 x
  │                                              ::::::::::          
  │  slot connections:                           ::::::::::          
  │  0¤ <─> ¤0                                   ::::::::::          
  │  1¤ <─> ¤1                                   ::::::::::          
  │  2¤ <─> ¤2                                   ::::::::::          
  │                                              ::::::::::          
  │ ├── Heun                                     ::::::::::          
  │ │  data slots:                               ::::::::::          
  │ │  [a] solution.membrane/V                   ├÷÷÷÷÷÷÷÷÷────── ¤0<─────┐
  │ │  [a] solution.sodium_channel_m_gate/m       ├÷÷÷÷÷÷÷÷m_gate ¤1<───┐ │
  │ │  [a] solution.sodium_channel_h_gate/h        ├÷÷÷÷÷÷÷h_gate ¤2<─┐ │ │
  │ │  [a] solution.potassium_channel_n_gate/n      ├÷÷÷÷÷÷────── ¤3 x│ │ │
  │ │  [a] additionalFieldVariable0                  ├÷÷÷÷÷─── aa ¤4 x│ │ │
  │ │  [a] additionalFieldVariable1                   ├÷÷÷÷─── bb ¤5 x│ │ │
  │ │  [a] leakage_current/i_L                         ├÷÷÷────── ¤6 x│ │ │
  │ │                                                   :::           │ │ │
  │ │ └── CellmlAdapter                                 :::           │ │ │
  │ └                                                   :::           │ │ │
  │                                                     :::           │ │ │
  │ ├── CrankNicolson                                   :::           │ │ │
  │ │  data slots:                                      :::           │ │ │
  │ │  [a] solution                                     ├÷÷─── vm ¤0<─┼─┼─┘
  │ │  [a] additionalFieldVariable0                      ├÷m_gate ¤1<─┼─┘
  │ │  [a] additionalFieldVariable1                       ├──── h ¤2<─┘
  │ │                                                                
  │ │ ├── FiniteElementMethod                                        
  │ │ │  data slots:                                                 
  │ │ │  [a] solution                                          vm ¤0 x
  │ │ │                                                              
  │ └                                                                
  └                                                                  
                                                                     
  Connection Types:
    +··+   Internal connection, no copy
    ════   Reuse variable, no copy
    ───>   Copy data in direction of arrow
    ─m──   Mapping between different meshes

  Referenced Meshes:
    [a] "MeshFiber", 1D regular fixed, linear Lagrange basis

With this example, the threre mechanisms to connect data slots can be seen:

* For connecting the variable ``solution.membrane/V`` in the `Heun` scheme to the ``solution`` variable in the `CrankNicolson` scheme, the connection was specified in the StrangSplitting scheme under the options ``"connectedSlotsTerm1To2"`` and ``"connectedSlotsTerm2To1"``. 

.. code-block:: python

    "connectedSlotsTerm1To2":     {0:0},   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
    "connectedSlotsTerm2To1":     {0:0},   # transfer the same back, in order to reuse field variables
    
Note that the slots do not have to have slot names defined.
    
* Connection of the `h` gating variable slots was done with the global setting "connectedSlots" as follows:

.. code-block:: python

  config = {
    ...
    "connectedSlots": [
      ("h", "h_gate"),      # connect the additional field variable in the output writer
      ("h_gate", "h"),
    ],
    ...
  
Here, the two slots ``h_gate`` and ``h`` are connected, ``h_gate`` is the name of the slot at the `CellmlAdapter` and ``h`` is the slot name at the additional field variable in the `CrankNicolson` scheme.

* The mechanism to connect slots by naming the slots the same is in use for the `m_gate` variable.


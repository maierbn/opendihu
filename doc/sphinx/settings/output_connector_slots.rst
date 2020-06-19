Output Connector Slots
===================================

Each solver has multiple *output connector slots* where field variables are exposed to surrounding solvers. 
When using a :doc:`/settings/coupling` or :doc:`/settings/splitting` scheme, data will be transferred between the two involved solvers.

After solution of term 1 is complete, all slots of term 2 will get the values of the connected slots from term 1. After the solution of term 2 is complete, the opposite connections will be used.

The connections between the slots of two solvers has to be specified in the settings.
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



Mesh
=======


inputMeshIsGlobal
^^^^^^^^^^^^^^^^^^
It specifies whether the given values and degrees of freedom are interpreted as local values or global values in the context of a parallel execution on multiple processes. It has no effect for serial execution.
It applies to all values given as mesh properties, such as node positions, element and node numbers, the physicalExtent, the number of elements, etc.

* If set to ``True``, all specified values and degrees of freedom are interpreted with global indexing. In this case, the same values should be given on all processes. Consequently, the program can be run on different numbers of processes with the same settings.
* If set to ``False``, all specified values and degrees of freedom are interpreted to be for the local portion of the own process, only.
  In parallel execution, each process has to get only its own range of values, which are typically different on each process. 

  For example, the number of elements is only specified for the local portion. Opendihu will compute the global number of elements from the local numbers on all the processes.

To provide different values for different MPI ranks, the own MPI rank number can be retrieved in the python settings. 
The last two command line arguments that are available in the python settings script are the own MPI rank number and the total number of ranks.

The advantage of the local specification is that each process only has to know its own portion of the whole problem. Internally there is no transfer of the local information to other processes. 
Thus, large problems can be computed with a high number of processes, where the global problem data would be too big to be stored by a single process.

The following example uses such a local specification. It sets the right hand side value of the last degree of freedom on the last MPI rank to 1.0 and all other values to 0.0.

.. code-block:: python

  # get own MPI rank number and number of MPI ranks
  rank_no = (int)(sys.argv[-2])
  n_ranks = (int)(sys.argv[-1])
  
  config = {
    "FiniteElementMethod" : {
      "inputMeshIsGlobal": False,
      "rightHandSide": {-1: 1.0} if rank_no == n_ranks-1 else {},
      
      # further options of FiniteElementMethod
      # ...
    }
  }

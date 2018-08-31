Mesh
  "nElements": n,
  "nodeDimension": 1,
  "nodePositions": [0, 1, ],  # as [x1,x2,x3,x4,...],  then set nodeDimension to 1/2/3. The number of entries must match the expected number of nodes (from nElements and mesh type)
  "nodePositions": [[0, 1], [2,3], [4,5,6], ...],
  "physicalExtent": [5.0, 8.0, 2.0]     # this can be used to set nodePositions
  "physicalExtent": 4.0,





"MultipleInstances": {
    "nInstances": nInstances,
    "instances": [
      {
        "ranks": [0,1,2],
        "GodunovSplitting": {
        }
      },
      {
      }
    ]
    "OutputWriter" : [],
  

Common names:
Element: has several nodes
Node: has several dofs, for Lagrange ansatz functions there is one dof per node, for cubic Hermite there are 2^D dofs per node (e.g. 1D: 2 dofs per node for value and derivative)
Dof: "Degree of freedom", an element of the function space that is spanned by the basis functions on the mesh. Dofs are located on the nodes, every node has at least 1 dof or more.
      To define a scalar field on the mesh one has to specify a value for each dof. For a vector field with nComponents components, (e.g. the geometryField with nComponents=3) there are
      nComponents*nDofs values needed.
Unknown: 
      A component of a dof. If only scalar fields are considered, nUnknowns = nDofs, otherwise nUnknows = nDofs*nComponents
      Note that the total number of unknowns for which the system has to be solved can be a multiple of the number of dofs when there are multiple components. 
      (E.g. 3D finite elasticity with 3-component displacement and 1-component pressure has 4 components per dofs, i.e. nUnkowns = nDofs*4)
      
No:   abbreviation for "number", in the sense of a counter value (i.e. not meaning "total sum" but as in numbering)
n:    E.g. nUnknowns, nDofs: "number", but this time meaning the total number/count of items

Consistent naming of local and global quantities:
<x>LocalWithGhosts
<x>LocalWithoutGhosts
<x>Global
E.g. nNodesLocalWithGhosts() or nDofsGlobal()

Numbering Schemes:
1. Global natural numbering, starting at 0 at the front left bottom corner, then continuing in x-direction, then y-direction, then z-direction (if in 3D)
   This numbering is used for global input of Dirichlet boundary conditions or right hand side values
2. Local numbering, per partition, starting from 0 to nDofsLocalWithoutGhosts-1, only numbering the non-ghost local dofs, in the same order as the global natural numbering but in the local partition.
   Then from nDofsLocalWithoutGhosts to nDofsLocalWithGhosts-1 the ghost dofs are numbered, again in the order given by the global natural numering.
   This numbering is used for accessing local vectors and matrices. One can access all values including ghosts by using the whole range or only access the non-ghost values by using the range [0,nDofsLocalWithoutGhosts-1]
3. Global Petsc numbering, starts at 0 on rank 0 and follows the local numbering for the non-ghost dofs, then it continues on rank 1 and follows the non-ghost dofs there and so on.
   This is the numbering that has to be used for accessing global Petsc vectors and matrices, however this is not needed because one can access these vectors and matrices through the local vectors and matrices.
   It is needed for the creation of the Vecs and Mats.

Howto debug:
1.GDB
define the following alias:
  alias gdb='gdb -ex=run --args '
then simply run with gdb:
gdb ./executable ../settings.py <further-arguments>
After the program crashes you can inspect the stacktrace with the command "bt" inside gdb.

2.Debugging output
- run with "-v" to enable all verbose output
- run with --v=1 or --v=2 etc. to enable verbose output to a given level
- run with -vmodule=partitioned_petsc_vec_structured.tpp=2,01_mesh_partition_structured.tpp=1 to enable verbose output of level2 only in the file partitioned_petsc_vec_structured.tpp and verbose output level 1 only in file 01_mesh_partition_structured.tpp
  Also wildcards (*) can be used, e.g. -vmodule=*vec*=3,*mat*=5,*mesh_partition*=1, then all files matching *vec*, *mat* or *mesh_partition* will get the specified output verbosity
  
3.Debugging parallel programs
- run program with mpirun and with "-pause" argument, example:
  mpirun -n 2 ./2_ranks_tests -v --gtest_filter=LaplaceTest.Structured1DHermite -pause
  Then it will stop with the following message:
0/2 INFO : Rank 0, PID 16614 is waiting for gdbResume=0 to become 1 

sudo gdb -p 16614

select-frame 2
set var gdbResume = 1
info locals 
continue
1/2 INFO : Rank 1, PID 16615 is waiting for gdbResume=0 to become 1 

sudo gdb -p 16615

select-frame 2
set var gdbResume = 1
info locals 
continue

  - now in two separate shell windows, execute "sudo gdb -p 16614" and "sudo gdb -p 16615". This attaches gdb to the two mpi processes. Inside gdb run the displayed commands "select-frame 2", "set var gdbResume = 1", etc. a
    After "continue" in both attached shells the program will continue. When it crashes, use "bt" to inspect the location again.
    
4.Memcheck
For segmentation faults that cannot be debugged with gdb, you can use valgrind with memcheck:
valgrind --tool=memcheck ./executable
There are a lot of "false positives" at the beginning while the python settings script is run. This is due to the python library overloading functions of memory management. Watch out for errors after these big outputs.

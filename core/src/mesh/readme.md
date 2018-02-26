This directory contains basic mesh data structures. There are 3 different mesh types:
* *structured regular fixed* which is a regular mesh with the same fixed mesh width for all coordinate directions. If it is used for finite element analysis the stiffness matrix will be set using precomputed stencils. This mesh type cannot be used for structural mechanics caculations with reference and actual configurations.
* *structured deformable* which is still regular in index space but can have arbitrary node positions. Contrary to the regular fixed mesh, this mesh can be deformed and used for structural mechanics.
* *unstructured deformable*. This is an unstructured cartesian grid.

Throughout the code not these basic classes are used but the BasisOnMesh classes which derive from one of these mesh classes and add functionality of a basis function to it.

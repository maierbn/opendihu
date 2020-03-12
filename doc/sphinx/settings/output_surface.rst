OutputSurface
===============

This class only modifies the output of OutputWriters. If there are 3D meshes present, it creates a 2D mesh on the surface of the 3D mesh and outputs only the 2D data.
This helps to reduce the output data size, usually what happens on the surface is the only interesting thing. Especially for problems with larger mesh sizes, the 3D mesh is often not practicable to store and visualize. Here, the 2D meshes can help.

C++ code:

.. code-block:: c

  OutputWriter::OutputSurface<
    // nested solver
  >

Python settings:

.. code-block:: python

  "OutputSurface": {
    "OutputWriter": [
      # output writer, e.g.:
      {"format": "Paraview", "outputInterval": (int)(1./dt_1D*output_timestep), "filename": "out/surface", "binary": True, "fixedFormat": False, "combineFiles": True},
    ],
    "face": ["1-"],
    
    # settings of the nested solver
  }
  
OutputWriter
--------------

The output writers that will output the 2D mesh only, not the 3D mesh.

face
--------
The 2D faces of the 3D volume which will be extracted. This has to be a list of faces, even if only a single face is required (then it is a list with one entry).
Possible values are `0-`, `0+`, `1-`, `1+`, `2-`, `2+`. :numref:`faces` shows the meaning for 1D, 2D and 3D meshes. In this case the 3D visualization (right picture) is relevant.o

.. _faces:
.. figure:: images/faces.svg
  :width: 100%

.. note::
  If the code runs in parallel and composite meshes are used, only the last sub mesh of the composite mesh is considered. Then each subdomain of that submesh has to span all us
  The reason is that each created 2D surface mesh has to take part in the collective communication functions of MPI. This is not the case if some surface meshes only need a subset of the MPI ranks of the other surface meshes.
  
  If unsure, try out the code. If errors or a crash occurs, then don't use `OutputSurface` in this context. 
  
  (This is the reason why the ``static_biceps_emg`` example only works in parallel either with only one face or with two faces but only with partitioning in z direction.)
  
  

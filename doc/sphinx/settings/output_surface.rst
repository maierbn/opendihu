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
    "face": "1-",
    
    # settings of the nested solver
  }
  
OutputWriter
--------------

The output writers that will output the 2D mesh only, not the 3D mesh.

face
--------
The 2D face of 3D volume which will be extracted. Possible values are `0-`, `0+`, `1-`, `1+`, `2-`, `2+`. :numref:`faces` shows the meaning for 1D, 2D and 3D meshes. In this case the 3D visualization (right picture) is relevant.o

.. _faces:
.. figure:: images/faces.svg
  :width: 100%

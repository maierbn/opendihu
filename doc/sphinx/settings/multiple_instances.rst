MultipleInstances
==================

.. code-block:: python

  "MultipleInstances": {
    "nInstances": nInstances,
    "instances": [
      {
        "ranks": [0,1,2],
        
        # further solver property objects, e.g. for GodunovSplitting
        "GodunovSplitting": {
        }
      },
      {
      ...
      }
    ]
    "OutputWriter" : [...],
  }
  

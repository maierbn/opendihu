# Script to be loaded in Paraview python console, ParaView 5.6.0-RC1
# It sets the coloring of every loaded mesh to "solution", useful for the fibre simulation where you would have to click 100 times on every fibre.
# Simply open the Python Shell inside paraview and execute
# >>> import set_fiber_coloring
#
# if the scripts needs to be reloaded (after a change), execute
# >>> reload(set_fiber_coloring)
#
# How to generate traces:
# import paraview.smstate as st
# st.smtrace.start_trace()
#    # do something
# st.smtrace.get_current_trace_output()   # this gets the current trace
# # write to file
# with open("file.py", "wb") as f:
#   f.write(st.smtrace.get_current_trace_output())
#
# # show the current working directory to know where the file was written
# import os
# print os.getcwd()
#

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# https://www.paraview.org/ParaView/Doc/Nightly/www/py-doc/paraview.simple.html

# loop over all sources
sources = GetSources()
for key,source in sources.iteritems():
    
  # set active source
  SetActiveSource(source)
  
  # get active view
  renderView1 = GetActiveViewOrCreate('RenderView')
  # uncomment following to set a specific view size
  # renderView1.ViewSize = [1477, 229]
  
  # get display properties
  sourceDisplay = GetDisplayProperties(source, view=renderView1)
  
  # change representation type
  sourceDisplay.SetRepresentationType('Wireframe')
  
  # Properties modified
  sourceDisplay.LineWidth = 2.0
  
  # set scalar coloring
  ColorBy(sourceDisplay, ('POINTS', 'solution'))
  
  # rescale color and/or opacity maps used to include current data range
  sourceDisplay.RescaleTransferFunctionToDataRange(True, False)
  
  # show color bar/color legend
  sourceDisplay.SetScalarBarVisibility(renderView1, True)

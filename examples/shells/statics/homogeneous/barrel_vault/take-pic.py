# trace generated using paraview version 5.13.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
ql = 'q'
nc = 2
quantity = ql + str(nc)
data = XMLUnstructuredGridReader(registrationName='data', FileName=['scolo_q4rs-hard-biased-512-' + ql + '.vtu'])

# Properties modified on data
data.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
dataDisplay = Show(data, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
dataDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

#changing interaction mode based on data extents
renderView1.InteractionMode = '3D'

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

#change interaction mode for render view
renderView1.InteractionMode = '2D'

# reset view to fit data bounds
renderView1.ResetCamera(0.0, 16.06969024216348, 0.0, 25.0, -5.84888892202555, 0.0, True, 0.9)

# set scalar coloring
ColorBy(dataDisplay, ('POINTS', quantity))

# rescale color and/or opacity maps used to include current data range
dataDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
dataDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for quantity
q1LUT = GetColorTransferFunction(quantity)

# get opacity transfer function/opacity map for quantity
q1PWF = GetOpacityTransferFunction(quantity)

# get 2D transfer function for quantity
q1TF2D = GetTransferFunction2D(quantity)

# get color legend/bar for q1LUT in view renderView1
q1LUTColorBar = GetScalarBar(q1LUT, renderView1)

# Properties modified on q1LUTColorBar
q1LUTColorBar.TitleFontSize = 44
q1LUTColorBar.LabelFontSize = 37

# reset view to fit data bounds
renderView1.ResetCamera(0.0, 16.06969024216348, 0.0, 25.0, -5.84888892202555, 0.0, True, 0.9)

# reset view to fit data bounds
renderView1.ResetCamera(0.0, 16.06969024216348, 0.0, 25.0, -5.84888892202555, 0.0, True, 0.9)

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(775, 777)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [2.713618148511602, 12.607138798105439, 55.730979483360834]
renderView1.CameraFocalPoint = [2.713618148511602, 12.607138798105439, -2.924444461012775]
renderView1.CameraParallelScale = 13.87447435465438

# save screenshot
SaveScreenshot(filename='snap' + quantity + '.png', viewOrLayout=renderView1, location=16, ImageResolution=[775, 777])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(775, 777)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [2.713618148511602, 12.607138798105439, 55.730979483360834]
renderView1.CameraFocalPoint = [2.713618148511602, 12.607138798105439, -2.924444461012775]
renderView1.CameraParallelScale = 13.87447435465438


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
# trace generated using paraview version 5.5.1-4-g5bcfde9

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PLOT3D Reader'
gridg = PLOT3DReader(QFileName=[''],
    FileName='C:\\Users\\yato\\Desktop\\ueno_flow\\grid.g',
    FunctionFileName='')
gridg.Functions = [110]

# set active source
SetActiveSource(gridg)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1154, 810]

# show data in view
gridgDisplay = Show(gridg, renderView1)

# trace defaults for the display properties.
gridgDisplay.Representation = 'Outline'
gridgDisplay.ColorArrayName = ['POINTS', '']
gridgDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
gridgDisplay.SelectOrientationVectors = 'None'
gridgDisplay.ScaleFactor = 40.0
gridgDisplay.SelectScaleArray = 'None'
gridgDisplay.GlyphType = 'Arrow'
gridgDisplay.GlyphTableIndexArray = 'None'
gridgDisplay.GaussianRadius = 2.0
gridgDisplay.SetScaleArray = [None, '']
gridgDisplay.ScaleTransferFunction = 'PiecewiseFunction'
gridgDisplay.OpacityArray = [None, '']
gridgDisplay.OpacityTransferFunction = 'PiecewiseFunction'
gridgDisplay.DataAxesGrid = 'GridAxesRepresentation'
gridgDisplay.SelectionCellLabelFontFile = ''
gridgDisplay.SelectionPointLabelFontFile = ''
gridgDisplay.PolarAxes = 'PolarAxesRepresentation'
gridgDisplay.ScalarOpacityUnitDistance = 3.0017744029070044
gridgDisplay.SelectInputVectors = [None, '']
gridgDisplay.WriteLog = ''

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
gridgDisplay.DataAxesGrid.XTitleFontFile = ''
gridgDisplay.DataAxesGrid.YTitleFontFile = ''
gridgDisplay.DataAxesGrid.ZTitleFontFile = ''
gridgDisplay.DataAxesGrid.XLabelFontFile = ''
gridgDisplay.DataAxesGrid.YLabelFontFile = ''
gridgDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
gridgDisplay.PolarAxes.PolarAxisTitleFontFile = ''
gridgDisplay.PolarAxes.PolarAxisLabelFontFile = ''
gridgDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
gridgDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

# set active source
SetActiveSource(gridg)

# Properties modified on gridg
gridg.QFileName = ['C:\\Users\\yato\\Desktop\\ueno_flow\\phys.q']
gridg.FunctionFileName = ''

# Rescale transfer function
gridgDisplay.ScaleTransferFunction.RescaleTransferFunction(0.677035947289, 1.00008823429)

# Rescale transfer function
gridgDisplay.OpacityTransferFunction.RescaleTransferFunction(0.677035947289, 1.00008823429)

# show data in view
gridgDisplay = Show(gridg, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(gridgDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
gridgDisplay.SetScalarBarVisibility(renderView1, True)

# Properties modified on gridgDisplay
gridgDisplay.SetScaleArray = [None, 'Density']

# Properties modified on gridgDisplay
gridgDisplay.OpacityArray = [None, 'Density']

# Properties modified on gridgDisplay
gridgDisplay.OSPRayScaleArray = 'Density'

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# Properties modified on gridgDisplay
gridgDisplay.SelectInputVectors = ['POINTS', 'Momentum']

# Properties modified on gridg
gridg.Functions = [110, 111, 112, 113, 120, 130, 140, 144, 153, 170, 184, 200, 201, 210, 211, 212]

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(gridgDisplay, ('POINTS', 'Pressure'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(vtkBlockColorsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
gridgDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
gridgDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Pressure'
pressureLUT = GetColorTransferFunction('Pressure')

# change representation type
gridgDisplay.SetRepresentationType('Surface')

# get color legend/bar for pressureLUT in view renderView1
pressureLUTColorBar = GetScalarBar(pressureLUT, renderView1)

# Properties modified on pressureLUTColorBar
pressureLUTColorBar.WindowLocation = 'AnyLocation'

# change scalar bar placement
pressureLUTColorBar.Position = [0.8648180242634316, 0.3320987654320988]
pressureLUTColorBar.ScalarBarLength = 0.33000000000000007

# current camera placement for renderView1
renderView1.CameraPosition = [100.0, -4.464163780212402, 994.1153045103262]
renderView1.CameraFocalPoint = [100.0, -4.464163780212402, 212.42845916748047]
renderView1.CameraParallelScale = 202.31544288083714

# get layout
layout1 = GetLayout()

# save screenshot
SaveScreenshot('C:/Users/yato/Desktop/ueno_flow/pressure.png', layout1, ImageResolution=[1154, 810],
    TransparentBackground=1)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# Properties modified on pressureLUTColorBar
pressureLUTColorBar.TitleOpacity = 0.0
pressureLUTColorBar.LabelOpacity = 0.0

# Properties modified on pressureLUTColorBar
pressureLUTColorBar.DrawAnnotations = 0

# Properties modified on pressureLUTColorBar
pressureLUTColorBar.ScalarBarThickness = 22

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [100.0, -4.464163780212402, 994.1153045103262]
renderView1.CameraFocalPoint = [100.0, -4.464163780212402, 212.42845916748047]
renderView1.CameraParallelScale = 202.31544288083714

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
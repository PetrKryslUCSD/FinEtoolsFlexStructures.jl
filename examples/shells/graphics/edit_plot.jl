using FinEtools
using FinEtools.MeshQuadrilateralModule
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using VisualStructures: plot_nodes, plot_midsurface, render, plot_space_box, plot_solid, space_aspectratio, default_layout_3d, plot_from_json
using PlotlyJS
using JSON

L = 42.0
H = 21.0


fens, fes = T3block(L, H, 7, 5)
fens.xyz = xyz3(fens)
sfes = FESetShellT3()
accepttodelegate(fes, sfes)
plots = cat(plot_nodes(fens),
    plot_midsurface(fens, fes; facecolor = "rgb(15, 155, 225)", lwidth = 4),
    dims = 1)
layout = default_layout_3d(autosize=true)
# layout = default_layout_3d(width=400, height=400)
layout[:showLegend] = true

config  = PlotConfig(plotlyServerURL="https://chart-studio.plotly.com", showLink=true)
pl = plot(plots, layout; config = config)

display(pl)

true


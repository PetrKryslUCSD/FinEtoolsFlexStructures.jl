using FinEtools
using FinEtoolsFlexStructures.FESetCorotBeamModule: FESetL2CorotBeam
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionCircle, CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, default_layout_3d, plot_from_json
using PlotlyJS
using JSON

L = 42
xyz = [0 0 0;
0 L/4 L*1/4;
L/4 L/4 L*2/4;
L/4 0 L*3/4;
0 0 L]
nL = 20

cs = CrossSectionCircle(s -> 5.9910, s -> [0.0, 0.0, 1.0])
fens, fes = frame_member(xyz, nL, cs)
plots = cat(plot_nodes(fens),
    plot_midline(fens, fes; color = "rgb(155, 155, 255)", lwidth = 4),
    dims = 1)
layout = default_layout_3d(width = 500, height = 500)
layout[:title] = ""
config  = PlotConfig(plotlyServerURL="https://chart-studio.plotly.com", showLink=true)



pl = plot(plots, layout; config=config)

display(pl)

savejson(pl, "plot.json")

# The below works when run interactively from the command line:
pl = plot_from_json("plot.json")

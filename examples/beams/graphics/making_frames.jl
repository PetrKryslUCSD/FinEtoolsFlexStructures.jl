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


pl = plot(plots, layout; 
    config=PlotConfig(showLink = true, toImageButtonOptions = Dict(:format=>"webp")))
display(pl)

f(frame) = if frame < 10; "0$frame"; else "$frame"; end

sleep(0.5)
let
    frame = 1
    for eyex in 1.0:-0.1:-1.0
        pl.plot.layout["scene"][:camera][:eye][:x] = eyex
        pl.plot.layout[:title] = "$(eyex)"
        react!(pl, pl.plot.data, pl.plot.layout)
        sleep(0.05)
        savefig(pl, "a$(f(frame)).png")
        frame += 1
    end
end

# Convert to a gif
# magick a*.png a.gif



# The below works when run interactively from the command line:
# savejson(pl, "plot.json")

# pl = plot_from_json("plot.json")

# sleep(0.5)
# for eyex in 1.0:-0.1:-1.0
#     pl.plot.layout["scene"][:camera][:eye][:x] = eyex
#     pl.plot.layout[:title] = "$(eyex)"
#     react!(pl, pl.plot.data, pl.plot.layout)
#     sleep(0.5)
#     savefig(pl, "a.png")
# end


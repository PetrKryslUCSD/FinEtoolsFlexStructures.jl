# # GARTEUR SM-AG19 Testbed: Construction of the geometry

# Source code: [`garteur_geometry_vis_tut.jl`](garteur_geometry_vis_tut.jl)

# ## Description

# This virtual test application is based on the test article  
# used by the GARTEUR Structures & Materials Action Group 19  
# which organized a Round Robin exercise where 12 European laboratories  
# tested a single structure between 1995 and 1997. The benchmark structure   
# was a laboratory structure built to simulate the dynamic behaviour  
# of an aeroplane. The structure was initially built for a benchmark  
# study on experimental modal analysis conducted by the  
# Structures and Materials Action Group (SM-AG19) of the Group  
# for Aeronautical Research and Technology in EURope (GARTEUR).  
# The test-bed was designed and manufactured by ONERA, France.

# ![](IMAC97photo.png)

# ### References

# - [GARTEUR] Ground Vibration Test Techniques, compiled by A Gravelle, GARTEUR
#   Structures & Materials Action Group 19 Technical report TP-115, 1999.
# - [BW] Etienne Balmes, Jan R. Wright, GARTEUR GROUP ON GROUND VIBRATION
#   TESTING | RESULTS FROM THE TEST OF A SINGLE STRUCTURE BY 12 LABORATORIES IN
#   EUROPE, Proceedings of DETC'97, 1997 ASME Design Engineering Technical
#   Conferences, September 14-17, 1997, Sacramento, California.

# ## Goals

# - Show how to construct model from multiple connected beams.
# - Demonstrate the use of massless connectors.
# - Visualize the structure interactively.
# 

# ## Geometry of the testbed airplane.

# The aluminum testbed was a rather simple structure which was reasonably
# dynamically representative of a simple airplane structure [GARTEUR](@ref
# References). It was composed of several beams simulating a fuselage with wings
# and a tail. Wing tip drums allowed to adjust bending and torsion frequencies
# similarly to airplane ones, with some very close modal frequencies. 

# ![](garteur-geom.png)

# The script included below defines the geometry of the structure, the
# cross-sectional properties, the connectivity, and the location of the nodes.

include("garteur_geometry_tut.jl")

# ## Basic visualization

# Here we use the `PlotlyJS` plotting library, with some simple utilities for
# generating geometry of beams.
using PlotlyJS
using VisualStructures: plot_solid, plot_space_box, render, default_layout_3d, save_to_json

# The colors are used to help distinguish between the individual parts of the model.
colors = [
"rgb(125, 155, 155)",  # 1 body
"rgb(125, 155, 155)",  # 2 wing
"rgb(125, 155, 155)",  # 3 drums
"rgb(125, 155, 155)",  # 4 vertical tail
"rgb(125, 155, 155)",  # 5 horizontal tail
"rgb(125, 155, 155)",  # 6 constraining plate
"rgb(15, 15, 155)",  # 7 massless connectors: wing beam to plate
"rgb(125, 15, 15)",  # 8 massless connectors: wing beam to drums
"rgb(125, 175, 15)",  # 9 massless connectors: fuselage to wing beam
"rgb(125, 175, 15)",  # 10 massless connectors: fuselage to tail
"rgb(15, 155, 155)",  # 11 massless connectors: sensors, point masses
]

# The geometry is defined in terms of "traces" (a trace is in the `PlotlyJS`
# parlance a graphical object; it could be a curve, surface, and a lot of other
# things). It is the name for a graphical object in `PlotlyJS`. The two points
# below define a box, which is helpful when setting the extents of the graphics
# display.
tbox = plot_space_box([[-1.2 * L -1.2 * L -1.2 * L]; [+1.2 * L +1.2 * L +1.2 * L]])
# For each finite element set in the array `fesa`, generate the graphics to
# represent that object.
traces = let traces = tbox
    for fes in fesa
        labl  = fes.label[1]
        tm = plot_solid(fens, fes; facecolor=colors[labl]);
        traces = cat(traces, tm; dims = 1)
    end
    traces
end
# The layout of the plot is defined with simple defaults.
layout = default_layout_3d(;width=900, height=900)
# Next, the graphics is rendered, and may be interacted with by zooming,
# panning, etc.
pl = render(traces; layout = layout)

# ## Visualizing the nodes

# In order to be able to discern the nodes we will reduce the opacity of the
# surfaces representing the beams, otherwise the nodes would be hidden by these
# surfaces. The geometry is defined in the same way as above.

tbox = plot_space_box([[-1.2 * L -1.2 * L -1.2 * L]; [+1.2 * L +1.2 * L +1.2 * L]])
traces = let traces = tbox
    for fes in fesa
        labl  = fes.label[1]
        tm = plot_solid(fens, fes; facecolor=colors[labl], opacity = 0.3);
        traces = cat(traces, tm; dims = 1)
    end
    traces
end
# Next we add to the "traces" the graphics representing all the nodes in the
# model as bright red dots.
using VisualStructures: plot_nodes
traces = cat(traces, plot_nodes(fens; color = "rgb(255, 15, 5)"); dims = 1)

# Finally, the graphics is presented.
layout = default_layout_3d(;width=900, height=900)
pl = render(traces; layout = layout)

nothing

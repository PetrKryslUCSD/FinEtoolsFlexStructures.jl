using PGFPlotsX
using CSV
using FinEtools

objects = []

f = CSV.File("$(@__DIR__)/spherical_cap_expl_22-KE-1.00000e+00.csv") 
xdata, ydata = f["t"][1:7:end],  f["v"][1:7:end];
@pgf p = PGFPlotsX.Plot(
{
color = "black",
line_width  = 1.0,
},
Coordinates([v for v in  zip(xdata, ydata)])
)
push!(objects, p)
push!(objects, LegendEntry("1.0"))

f = CSV.File("$(@__DIR__)/spherical_cap_expl_22-KE-1.00000e-01.csv") 
xdata, ydata = f["t"][1:33:end],  f["v"][1:33:end];
@pgf p = PGFPlotsX.Plot(
{
only_marks,
        color = "black",
        line_width  = 1.0,
        mark = "x",
},
Coordinates([v for v in  zip(xdata, ydata)])
)
push!(objects, p)
push!(objects, LegendEntry("0.1"))

f = CSV.File("$(@__DIR__)/spherical_cap_expl_22-KE-1.00000e+01.csv") 
xdata, ydata = f["t"][17:33:end],  f["v"][17:33:end];
@pgf p = PGFPlotsX.Plot(
{
only_marks,
        color = "black",
        line_width  = 1.0,
        mark = "+",
},
Coordinates([v for v in  zip(xdata, ydata)])
)
push!(objects, p)
push!(objects, LegendEntry("10.0"))

@pgf ax = Axis(
{
xlabel = "Time [ms]",
ylabel = "Kinetic energy [lbf*in]",
# ymin = -1.1e-5,
# ymax = 1.1e-5,
xmin = 0.0,
xmax = 1.0,
xmode = "linear", 
ymode = "linear",
yminorgrids = "true",
grid = "both",
legend_style = {
at = Coordinate(0.5, 1.05),
anchor = "south",
legend_columns = -1
},
},
objects...
)

display(ax)
pgfsave("spherical_cap_expl_22-KE.pdf", ax)


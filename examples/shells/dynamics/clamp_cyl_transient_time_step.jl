using PGFPlotsX
using CSV
using FinEtools

results = [
(1.00000e-02, Any[2.52679e-07, 2.52672e-07, 2.52662e-07, 2.52651e-07, 2.52671e-07, 2.52651e-07, 1.90903e-07, 1.34989e-07])                                                               
(2.00000e-03, Any[2.52633e-07, 2.52675e-07, 2.52699e-07, 2.52659e-07, 2.41211e-07, 2.15746e-07, 1.52555e-07, 1.07873e-07])                                                               
(1.00000e-03, Any[2.24037e-07, 2.24064e-07, 2.24032e-07, 2.24037e-07, 2.06195e-07, 1.84427e-07, 1.30409e-07, 9.22134e-08])                                                               
(1.00000e-04, Any[1.93031e-07, 1.93028e-07, 1.93035e-07, 1.93027e-07, 1.84096e-07, 1.64660e-07, 1.16432e-07, 8.23301e-08])                                                               
(1.00000e-05, Any[1.92667e-07, 1.92665e-07, 1.92671e-07, 1.92664e-07, 1.83822e-07, 1.64415e-07, 1.16259e-07, 8.22077e-08])  
]

drilling_stiffness_scales = [0.01, 0.1, 1.0, 1.5, 2.0, 2.5, 5.0, 10.0]

objects = []
marks = ["x", "o", "square", "+", "diamond"]

for i in eachindex(results)
    xdata, ydata = drilling_stiffness_scales, results[i][2] ./ maximum(results[i][2]);
    @pgf p = PGFPlotsX.Plot(
    {
    color = "black",
    mark = "$(marks[i])",
    line_width  = 1.0,
    },
    Coordinates([v for v in  zip(xdata, ydata)])
    )
    push!(objects, p)
    push!(objects, LegendEntry("t/R=$(results[i][1]/0.1)"))
end

@pgf ax = Axis(
{
xlabel = "Drilling Stiffness Multiplier [ND]",
ylabel = "Normalized time step [ND]",
# ymin = -1.1e-5,
# ymax = 1.1e-5,
# xmin = 0.0,
# xmax = 1.0,
xmode = "log", 
ymode = "linear",
yminorgrids = "true",
grid = "both",
legend_style = {
at = Coordinate(0.35, 0.05),
anchor = "south",
# legend_columns = -1
},
},
objects...
)

display(ax)
pgfsave("clamped_cylinder_vibration-time_step.pdf", ax)


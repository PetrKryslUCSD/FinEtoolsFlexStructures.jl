[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl/actions)
[![codecov](https://codecov.io/gh/PetrKryslUCSD/FinEtoolsFlexStructures.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PetrKryslUCSD/FinEtoolsFlexStructures.jl)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtoolsFlexStructures.jl/dev)


# FinEtoolsFlexStructures.jl

[`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl) used for 
- Simulations of large-displacement response of three-dimensional flexible-beam structures. Linear static analysis, modal analysis, linear buckling analysis. Nonlinear statics and dynamics;
- Simulations of shell structures. Linear static analysis, modal analysis.

## Current limitations

- Only elastic structures can be modeled.
- Only simple solid beam cross sections are implemented. Open thin-walled beams cannot be simulated at this point since the warping of the section is not enabled.
- Shells need to be homogeneous. The layered model suitable for laminates and such is not implemented yet.
- Neither beams nor shells can be attached with an offset (eccentricity).

## News

- 12/14/2021: The shell element T3FF fully tested, and described in a paper (to be submitted).
- 08/18/2021: Implemented fully linear statics and dynamics of shells using the DSG triangle with various improvements. Robust and accurate element.
- 05/23/2021: Updated for Julia 1.6.
- 08/23/2020: Added a separate tutorial package, [FinEtoolsFlexStructuresTutorials.jl](https://petrkryslucsd.github.io/FinEtoolsFlexStructuresTutorials.jl)).
- 08/16/2020: Described the tutorials.
- 07/27/2020: Added documentation.
- 02/25/2020: Beams: Nonlinear static analysis implemented.
- 02/20/2020: Beams: Nonlinear transient dynamic analysis implemented.
- 02/16/2020: Beams: Buckling analysis implemented.

## Examples

There are a number of examples, which may be executed as described in the conceptual guide to [`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl).

## Visualization

This is possible using [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) with the package [`FinEtoolsFlexStructures.VisUtilModule`](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.VisUtilModule.jl).

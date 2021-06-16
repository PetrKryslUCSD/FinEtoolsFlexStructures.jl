[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl.svg?branch=master)](https://travis-ci.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtoolsFlexStructures.jl/dev)

# FinEtoolsFlexStructures.jl

FinEtools used for the simulation of large-displacement response of three-dimensional flexible-beam structures. Linear static analysis, modal analysis, linear buckling analysis. Nonlinear statics and dynamics.

![](http://hogwarts.ucsd.edu/~pkrysl/site.images/circle-twist-anim.gif)

## News

- 05/23/2021: Updated for Julia 1.6.
- 08/23/2020: Added a separate tutorial package, [FinEtoolsFlexStructuresTutorials.jl](https://petrkryslucsd.github.io/FinEtoolsFlexStructuresTutorials.jl)).
- 08/16/2020: Described tutorials.
- 07/27/2020: Added documentation.
- 02/25/2020: Nonlinear static analysis implemented.
- 02/20/2020: Nonlinear transient dynamic analysis implemented.
- 02/16/2020: Buckling analysis implemented.

## Tutorials

Clone the package to your working directory:
```
git clone https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl.git
```

Change your working directory to `FinEtoolsFlexStructures`. Start Julia and run
the following:

```
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

The tutorials are in the form of Julia scripts in the folder `FinEtoolsFlexStructures.jl/docs/src/tutorials`. The markdown generated from these files is also in the same folder.

To view the markdown, follow the link to the documentation.

To run a tutorial, head over to the `tutorials` folder, open the tutorial script, and evaluate in Julia.

## Examples

There are a number of examples, which may be executed as described in the conceptual guide to [`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl).

## Visualization

This is possible with the package [`FinEtoolsFlexStructures.VisUtilModule`](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.VisUtilModule.jl).

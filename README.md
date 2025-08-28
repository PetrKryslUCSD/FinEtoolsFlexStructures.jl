[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl/actions)
[![codecov](https://codecov.io/gh/PetrKryslUCSD/FinEtoolsFlexStructures.jl/branch/main/graph/badge.svg?token=5MHDMHEFCY)](https://codecov.io/gh/PetrKryslUCSD/FinEtoolsFlexStructures.jl)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtoolsFlexStructures.jl/dev)


# FinEtoolsFlexStructures.jl

[`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl) used for 
- Simulations of large-displacement response of three-dimensional flexible-beam
  structures. Linear static analysis, modal analysis, linear buckling analysis.
  Nonlinear statics and dynamics;
- Simulations of large displacement response of three dimensional truss structures.
  Nonlinear statics, including path tracing.
- Simulations of shell structures. Linear static analysis, modal analysis,
  explicit dynamic analysis. Shells can be homogeneous or layered
  (laminated, composite).

## Current limitations

- Only elastic structures can be modeled.
- Only simple solid or closed thin-walled beam cross sections are implemented.
  Open thin-walled beams cannot be simulated at this point since the warping of
  the section is not enabled.
- Shell reference surface offset not yet tested.
- Neither general beams nor shells can be attached with an offset
  (eccentricity). There is a linear beam model that can be used with
  eccentricities (offsets of the attachment points relative to the nodes).

## News

- 08/28/2025: Update examples.


[Past news](#past-news)


## Tutorials

There are a number of tutorials explaining the use of this package.
Check out the [index](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl/blob/main/tutorials/index.md). The tutorials  can be executed as
follows:

- Download the package or clone it.
```
git clone https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl.git
```
- Change into the `tutorials` folder: `cd .\FinEtoolsFlexStructures.jl\tutorials`.
- Start Julia: `julia`.
- Activate the environment:
```
using Pkg; Pkg.activate("."); Pkg.instantiate();
```
- Execute the desired tutorial. Here we arbitrarily pick one:
```
include("circle_modal_tut.jl")
```

When running a tutorial with graphical plotting in VSCode, make sure to turn off
(uncheck) "Julia: Use Plot Pane" in the  Preferences. Otherwise the animations
will not work.


## Examples

Let us assume that the working directory   is `"FinEtoolsFlexStructures.jl"`, perhaps
as a result of cloning the repository.

The project needed for running examples can be activated and initialized by
```
using Pkg; Pkg.activate("./examples"); Pkg.instantiate();
```

There are a number of examples, which may be executed as described in the
conceptual guide to [`FinEtools`]
(https://github.com/PetrKryslUCSD/FinEtools.jl). As an example:

```
julia> include(".\\examples\\shells\\dynamics\\dcbs_vibration_examples.jl")                 
WARNING: replacing module dbcs_vibration_examples.             
[ Info: All examples may be executed with                      
using .Main.dbcs_vibration_examples; Main.dbcs_vibration_examples.allrun()       
                  
julia> using .Main.dbcs_vibration_examples; Main.dbcs_vibration_examples.allrun()                         
                  
#####################################################          
# test_convergence
[ Info: FV12 free vibration, formulation=FinEtoolsFlexStructures.FEMMShellT3FFModule                        
[ Info: Mesh: 1 elements per side
count(fens) = 143 
fs = [21.613635088738015, 29.23401281588026, 30.925491823018966, 34.36778956332975]                         
[ Info: Mesh: 2 elements per side   
count(fens) = 567 
fs = [20.569847475978634, 26.277349178877216, 30.027181006351178, 30.68607366768112]                        
[ Info: Mesh: 3 elements per side   
count(fens) = 2255
fs = [20.350441226155656, 25.67267626791537, 29.124064583761363, 30.620439988286456]                        
[ Info: Mesh: 4 elements per side   
count(fens) = 8991
fs = [20.301524870325565, 25.533290848730623, 28.914284995255777, 30.620822302876647]                       
true      
```
There is usually some indication of what the correct answer should be in 
the document string of the module: refer to the Julia file defining the examples.

## Visualization

This is possible with [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) with the package [`VisualStructures`](https://github.com/PetrKryslUCSD/VisualStructures.jl).
Static plots or animation of deformation during a static or dynamic simulation can be done.

Export to  [Paraview](https://www.paraview.org/)  is also available. Static
pictures and time collections (useful for animations) are supported this way.

## Documentation

- Paper on explicit dynamics accepted for publication  in the
  International Journal for Numerical Methods in Engineering. [Draft is
  available in PDF.](docs/expl-shells-compressed.pdf)
- [Paper](https://doi.org/10.1002/nme.6944) describing the robust
  triangular flat-facet shell element has been accepted for publication in the
  International Journal for Numerical Methods in Engineering. [Draft is
  available in PDF.](docs/shells-submitted.pdf)

## <a name="past-news"></a>Past news

- 01/06/2024: Update tutorials.
- 01/04/2024: Add a Riks path-tracing solver.
- 01/01/2024: Add a nonlinear corotational truss element. 
- 12/22/2024: Add vibration examples for angle-ply laminated composites.
- 12/19/2024: Add benchmark examples for laminated composites.
- 12/19/2024: Fix transformation bug in composite implementation. 
- 12/15/2024: Reorganize examples into separate homogeneous and layered shell folders.
- 12/10/2024: Update to Julia 1.11. Reorganize and update examples.
- 05/27/2024: Fix documentation of test module.
- 02/25/2024: Update for FinEtools 8.
- 12/30/2023: Update for Julia 1.10.
- 12/21/2023: Tutorials merged back into the package tree.
- 09/26/2023: Multiple examples revised. Numerous documentation improvements.
- 08/14/2023: Update for Julia 1.9.2. Revised nomenclature of layup transformation matrices.
- 07/13/2023: Update for Julia 1.9 and FinEtools 7.
- 05/26/2022: Paper on explicit dynamics accepted for publication  in the
  International Journal for Numerical Methods in Engineering. [Draft is
  available in PDF.](docs/expl-shells-compressed.pdf)
- 05/07/2022: Upgrade to Julia 1.7.2.
- 02/13/2022: [Paper](https://doi.org/10.1002/nme.6944) describing the robust
  triangular flat-facet shell element has been accepted for publication in the
  International Journal for Numerical Methods in Engineering. [Draft is
  available in PDF.](docs/shells-submitted.pdf)
- 01/29/2022: Explicit dynamics with CSR sparse matrix parallel multiplication.
- 12/31/2021: Implement model for layered (laminated, composite) plates and shells. 
- 12/20/2021: Reorganize examples into a project.
- 12/14/2021: The shell element T3FF fully tested, and described in a paper (submitted).
- 08/18/2021: Implemented linear statics and dynamics of shells using the DSG triangle with various improvements. Robust and accurate element.
- 05/23/2021: Update for Julia 1.6.
- 08/23/2020: Add a separate tutorial package, [FinEtoolsFlexStructuresTutorials.jl](https://petrkryslucsd.github.io/FinEtoolsFlexStructuresTutorials.jl)).
- 08/16/2020: Describe the tutorials.
- 07/27/2020: Add documentation.
- 02/25/2020: Beams: Nonlinear static analysis implemented.
- 02/20/2020: Beams: Nonlinear transient dynamic analysis implemented.
- 02/16/2020: Beams: Buckling analysis implemented.

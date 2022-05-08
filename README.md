[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl/actions)
[![codecov](https://codecov.io/gh/PetrKryslUCSD/FinEtoolsFlexStructures.jl/branch/main/graph/badge.svg?token=5MHDMHEFCY)](https://codecov.io/gh/PetrKryslUCSD/FinEtoolsFlexStructures.jl)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtoolsFlexStructures.jl/dev)


# FinEtoolsFlexStructures.jl

[`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl) used for 
- Simulations of large-displacement response of three-dimensional flexible-beam structures. Linear static analysis, modal analysis, linear buckling analysis. Nonlinear statics and dynamics;
- Simulations of shell structures. Linear static analysis, modal analysis, explicit dynamic analysis. Shells can be homogeneous or layered (laminated, composite).

## Current limitations

- Only elastic structures can be modeled.
- Only simple solid beam cross sections are implemented. Open thin-walled beams cannot be simulated at this point since the warping of the section is not enabled.
- Shell reference surface offset not yet tested.
- Neither beams nor shells can be attached with an offset (eccentricity).

## News

- 05/07/2022: Upgraded to Julia 1.7.2.
- 02/13/2022: [Paper](https://doi.org/10.1002/nme.6944) describing the robust triangular flat-facet shell element has been accepted for publication in the International Journal for Numerical Methods in Engineering. [Draft is available in PDF.](docs/shells-submitted.pdf)
- 01/29/2022: Explicit dynamics with CSR sparse matrix parallel multiplication.
- 12/31/2021: Implemented model for layered (laminated, composite) plates and shells. 
- 12/20/2021: Reorganized examples into a project.
- 12/14/2021: The shell element T3FF fully tested, and described in a paper (submitted).
- 08/18/2021: Implemented linear statics and dynamics of shells using the DSG triangle with various improvements. Robust and accurate element.
- 05/23/2021: Updated for Julia 1.6.
- 08/23/2020: Added a separate tutorial package, [FinEtoolsFlexStructuresTutorials.jl](https://petrkryslucsd.github.io/FinEtoolsFlexStructuresTutorials.jl)).
- 08/16/2020: Described the tutorials.
- 07/27/2020: Added documentation.
- 02/25/2020: Beams: Nonlinear static analysis implemented.
- 02/20/2020: Beams: Nonlinear transient dynamic analysis implemented.
- 02/16/2020: Beams: Buckling analysis implemented.

## Testing

The package `FinEtoolsFlexStructures.jl` is registered. Simply do
```
using Pkg; Pkg.add("FinEtoolsFlexStructures"); 
Pkg.test("FinEtoolsFlexStructures"); 
```

The package can also be cloned.
Let us assume that the working directory   is `"FinEtoolsFlexStructures.jl"`, 
which is the result of cloning the repository.
The present package can be tested with
```
using Pkg; Pkg.activate("."); Pkg.instantiate(); 
using Pkg; Pkg.test(); 
```

## Examples

The project was developed with Julia 1.6.4, and is operational with any version above that.

Let us assume that the working directory   is `"FinEtoolsFlexStructures.jl"`, perhaps
as a result of cloning the repository.

The project needed for running examples can be activated and initialized by
```
using Pkg; Pkg.activate("./examples"); Pkg.instantiate(); using Revise
```

There are a number of examples, which may be executed as described in the conceptual guide to [`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl). As an example:
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


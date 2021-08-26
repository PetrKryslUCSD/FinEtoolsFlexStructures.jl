using Documenter, FinEtoolsFlexStructures

makedocs(
	modules = [FinEtoolsFlexStructures],
	doctest = false, clean = true,
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "FinEtoolsFlexStructures.jl",
	pages = Any[
			"Home" => "index.md",
			"How to guide" => "guide/guide.md",
			"Reference" => "man/reference.md"	
		],
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl.git",
    devbranch = "main"
)

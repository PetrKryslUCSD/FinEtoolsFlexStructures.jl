using Documenter, FinEtoolsFlexStructures

makedocs(
	modules = [FinEtoolsFlexStructures],
	doctest = false, clean = true,
	warnonly = Documenter.except(:linkcheck, :footnote),
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "FinEtoolsFlexStructures.jl",
	pages = Any[
			"Home" => "index.md",
			"How to guide" => "guide/guide.md",
			"Reference" => "man/man.md"	
		],
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl.git",
    devbranch = "main"
)

using Literate

# Copy tutorial sources
for t in readdir("$(@__DIR__)")
    if occursin(r".*_tut.jl", t)
        println("\nTutorial $t in $(@__DIR__)\n")
        Literate.markdown(joinpath(@__DIR__, t), "$(@__DIR__)"; documenter=false);
        # cp(t, "../../../src/" * t, force = true)
    end
end

# Copy ancillary files
# for a in ["fast_top_ref.txt"]
#     println("\nAncillary $a in $(pwd())\n")
#     cp(a, "../../../src/" * a, force = true)
# end

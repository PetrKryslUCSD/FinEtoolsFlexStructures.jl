using Test
@time @testset "Beam section" begin
    include("test_section.jl")
end
@time @testset "Beam mesh" begin
    include("test_mesh.jl")
end
@time @testset "Beam modal" begin
    include("test_modal.jl")
end
@time @testset "Beam buckling" begin
    include("test_buckling.jl")
end
@time @testset "Beam transient dynamics" begin
    include("test_dyn.jl")
end
@time @testset "Shell statics" begin
    include("test_shell_statics.jl")
end
@time @testset "Shell dynamics" begin
    include("test_shell_dynamics.jl")
end

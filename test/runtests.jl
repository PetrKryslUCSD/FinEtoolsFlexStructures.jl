using Test

@time @testset "Composite shell statics Q4RS" begin
    include("test_composite_shell_statics_q4rs.jl")
end
@time @testset "Composite Layup" begin
    include("test_composite_layup.jl")
end
@time @testset "Composite shell dynamics T3FF" begin
    include("test_composite_shell_dynamics.jl")
end
@time @testset "Composite shell statics T3FF" begin
    include("test_composite_shell_statics.jl")
end
@time @testset "Shell statics Q4RS" begin
    include("test_q4rs_shell_statics.jl")
end
@time @testset "Beam transient dynamics" begin
    include("test_beam_dyn.jl")
end
@time @testset "Utilities" begin
    include("test_utilities.jl")
end
@time @testset "Beam section" begin
    include("test_section.jl")
end
@time @testset "Beam mesh" begin
    include("test_beam_mesh.jl")
end
@time @testset "Beam statics" begin
    include("test_beam_linear_statics.jl")
end
@time @testset "Beam modal" begin
    include("test_beam_modal.jl")
end
@time @testset "Beam buckling" begin
    include("test_beam_buckling.jl")
end
@time @testset "Shell dynamics" begin
    include("test_shell_dynamics.jl")
end

@time @testset "Shell resultants" begin
    include("test_shell_resultants.jl")
end
@time @testset "Shell statics" begin
    include("test_shell_statics.jl")
end

@time @testset "Rigid link" begin
    include("test_rigid.jl")
end
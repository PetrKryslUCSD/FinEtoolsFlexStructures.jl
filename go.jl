using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

# include(".\\examples\\nonlin_statics\\twisting_circle_examples.jl");                                                                         
# pl = twisting_circle_examples.twisting_circle_1() 
# savejson(pl, "plot.json")
# open("plot.json", "r") do f; plot(JSON.parse(Plot, String(read(f)))) end

# include(".\\examples\\nonlin_statics\\twisting_circle_examples.jl");                                                                         
# twisting_circle_examples.twisting_circle() 

# snapnum(n) = begin
#     s = ""
#     if n <= 9
#         s = "00$(n)"
#     elseif n <= 99
#         s = "0$(n)"
#     else
#         s = "$(n)"
#     end
#     s
# end
# run(`cp ../twi*.png .`)
# K = 1
# for i in 1:5:501
#     if mod(i, 10) == 1   
 #         global K += 1    
#     end
# end
# run(`magick snap*.png anim.gif`)

# include(".\\examples\\nonlin_transient\\twisting_circle_examples.jl");    
# twisting_circle_examples.twisting_circle() 

include(".\\examples\\nonlin_statics\\curved_cantilever_examples.jl");      
curved_cantilever_examples.curved_cantilever_thin() 

include(".\\examples\\meshing\\mesh_examples.jl");
mesh_examples.allrun()
mesh_examples.line_mesh_solid()

include(".\\examples\\modal\\beam_modal_examples.jl");
beam_modal_examples.allrun()

include(".\\examples\\modal\\argyr_l_frame_modal_examples.jl");
# argyr_l_frame_modal_examples.allrun()
argyr_l_frame_modal_examples.argyr_l_frame_modal_anim()

include(".\\examples\\linear_buckling\\tippling_examples.jl");
tippling_examples.allrun()

include(".\\examples\\nonlin_transient\\fast_top_examples.jl");
@time fast_top_examples.fasttop2()

include(".\\examples\\nonlin_transient\\slow_top_examples.jl");
@time slow_top_examples.slowtop3()

include(".\\examples\\nonlin_transient\\argyr_swing_examples.jl");          
argyr_swing_examples.argyr_swing_animated() 

include(".\\examples\\nonlin_transient\\argyr_swing_examples.jl");          
argyr_swing_examples.argyr_swing_compare() 

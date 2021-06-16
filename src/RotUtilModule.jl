module RotUtilModule

using LinearAlgebra: norm, mul!
using FinEtools

"""
    initial_Rfield(fens)

Set up initial rotation field (rotation matrices are all identities)
"""
function initial_Rfield(fens)
    f = NodalField(zeros(size(fens.xyz,1), 9))
    for i in 1:count(fens)
        f.values[i, 1] = f.values[i, 5] = f.values[i, 9] = 1.0
    end
    return f
end

"""
    update_rotation_field!(Rfield, dchi)

Update rotation field by exponential incremental rotation. 
"""
function update_rotation_field!(Rfield, dchi)
    Rs = Rfield.values; 
    Adelta = dchi.values[:, 4:6]; 
    R = fill(0.0, 3, 3)
    Rdelta = fill(0.0, 3, 3)
    Rupdated = fill(0.0, 3, 3)
    for i in 1:size(Rs,1)
        R[:] .= @view Rs[i, :];
        rotmat3!(Rdelta, @view Adelta[i,:])
        mul!(Rupdated, Rdelta, R)
        Rs[i,:] .= @view Rupdated[:];
    end 
    return Rfield
end

"""
    linear_update_rotation_field!(Rfield, dchi)

Update rotation field by linear incremental rotation. 
"""
function linear_update_rotation_field!(Rfield, dchi)
    Rs = Rfield.values; 
    Adelta = dchi.values[:, 4:6]; 
    R = fill(0.0, 3, 3)
    Sdelta = fill(0.0, 3, 3)
    Rupdated = fill(0.0, 3, 3)
    for i in 1:size(Rs,1)
        R[:] .= @view Rs[i, :];
        skewmat!(Sdelta, @view Adelta[i,:])
        mul!(Rupdated, Sdelta, R)
        Rs[i,:] .= @view Rupdated[:];
    end 
    return Rfield
end

end # module

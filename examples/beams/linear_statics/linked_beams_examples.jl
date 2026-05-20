"""
pystran - Python package for structural analysis with trusses and beams

(C) 2025, Petr Krysl, pkrysl@ucsd.edu

# Bending of a two incompressible cantilevers rigidly linked at the free end

## Problem description:

Two parallel cantilevers are rigidly linked at the free end.
They bend in sync due to the rigid link.

Displacements and bending moments are provided in the reference.

## References

This is the AFNOR SSLL05/89 test case.

Original source: "Guide de validation des progiciels de calcul de structures"
publié par l'AFNOR 1990 (ISBN 2-12-486611-7).

Data taken from: ICAB Force Exemples Exemples de calculs de statique pour ICAB
Force. www.icab.fr
"""
module linked_beams_examples

using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexStructures.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexStructures.FEMMLinBeamModule
using FinEtoolsFlexStructures.FEMMLinBeamModule: FEMMLinBeam
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtoolsFlexStructures.FEMMRigidLinkModule
using FinEtoolsFlexStructures.FEMMRigidLinkModule: FEMMRigidLink
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio, save_to_json, plot_local_frames

# The geometry of the problem is defined by the length of the beams and the
# distance between them.
const E = 200e9
# Original: ridiculously large cross section
# const A = 1.0
# const I = 4 / 3 * 10^-8
# For the cross section properties below, the python code gives these deflections:
# [ 4.62677504e-06 -1.52741910e-03 -2.77606502e-04]
# [-4.62677504e-06 -1.53241835e-03 -2.82230944e-04]
# Analytical bending moment magnitude:  500.0
const A = 0.01
const I = 4 / 3 * 10^-6
const H = sqrt(12 * I / A)
const B = A / H
const h = 2.0
const Q = 1.0e3 # applied force 
const d = 0.2


# The rigid link formulation requires a penalty for each degree of freedom.
# This penalty needs to be of appropriate magnitude, commensurate with the
# other stiffness in the system. This particular formulation does not require
# the penalty to be many times larger, just in the ballpark. Here we take a
# fraction of the `E*A`  product.
const Gamma = 1e8

function test(nel=1, visualize=true)
    @show H, B
    cs = CrossSectionRectangle(s -> B, s -> H, s -> [0.0, 0.0, 1.0]) # Bernoulli
    members = []
    xyz = [[0 0 0]; [h 0 0]]
    push!(members, frame_member(xyz, nel, cs))
    xyz = [[0 -d 0]; [h -d 0]]
    push!(members, frame_member(xyz, nel, cs))
    fens, fes = merge_members(members)

    material = MatDeforElastIso(DeforModelRed3D, E, 0.0)
    
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    clamped = selectnode(fens; box = [0 0 -d 0 0 0], inflate = d/1000)
    for i in [1, 2, 3, 4, 5, 6] # clamped
        setebc!(dchi, clamped, true, i)
    end
    for k in eachindex(fens)
        for i in [3, 4, 5] # constrained to XY plane
            setebc!(dchi, [k], true, i)
        end
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMLinBeam(IntegDomain(fes, GaussRule(1, 1), B * H), material)

    K = FEMMLinBeamModule.stiffness(femm, geom0, u0, Rfield0, dchi);
    # @show K_ff = Matrix(matrix_blocked(K, nfreedofs(dchi))[:ff])
    linked = selectnode(fens; box = [h h -d 0 0 0], inflate = h/1000)
    rb = FEMMRigidLink(
        IntegDomain(FESetP1(reshape(linked, length(linked), 1)), PointRule()), 
        linked[1], 
        Gamma)

    Krb = FEMMRigidLinkModule.stiffness(rb, geom0, u0, Rfield0, dchi)

    loaded = selectnode(fens; box = [h h -d -d 0 0], inflate = h/1000)
    loadbdry = FESetP1(reshape(loaded, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = fill(0.0, 6); q[2] = -Q
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve_blocked!(dchi, K + Krb, F)

    @show dchi.values[linked, :]

    if visualize
        scaling = 1e3
        dchi.values .*= scaling
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[-h -h -h];
         [h h h]]),
                    plot_nodes(fens),
                    plot_solid(fens, fes;
                               x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
                    dims = 1)
        pl = render(plots)
    end

    # # Using the stiffness coefficients of a beam in the basic configuration, we
    # # can calculate the deflection from the relationship between the shear force
    # # in the beams (there are two!) and the applied force:

    # w2 = F / (12 * E * I / h**3) / 2
    # print("Analytical deflection in the direction of the force: ", w2)
    # print("Expected tip deflection: ", [0.0, -1.25250329e-01, 0])
    # for j in m["joints"].values():
    #     print(j["displacements"])

    # plots.setup(m)
    # plots.plot_members(m)
    # plots.plot_deformations(m)
    # plots.show(m)

    # # Again using the basic stiffness, we deduce the correct magnitude of the bending moment:
    # M = 6 * h * E * I / h**3 * w2
    # print("Analytical bending moment magnitude: ", M)
    # plots.setup(m)
    # plots.plot_members(m)
    # plots.plot_bending_moments(m)
    # plots.show(m)
end # function test

function test_XZ(nel=1, visualize=true)
    @show H, B
    cs = CrossSectionRectangle(s -> B, s -> H, s -> [0.0, 1.0, 0.0]) # Bernoulli
    members = []
    xyz = [[0 0 0]; [h 0 0]]
    push!(members, frame_member(xyz, nel, cs))
    xyz = [[0 -d 0]; [h 0 -d]]
    push!(members, frame_member(xyz, nel, cs))
    fens, fes = merge_members(members)

    material = MatDeforElastIso(DeforModelRed3D, E, 0.0)
    
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    clamped = selectnode(fens; box = [0 0 0 0 -d 0], inflate = d/1000)
    for i in [1, 2, 3, 4, 5, 6] # clamped
        setebc!(dchi, clamped, true, i)
    end
    # for k in eachindex(fens)
    #     for i in [3, 4, 5] # constrained to XY plane
    #         setebc!(dchi, [k], true, i)
    #     end
    # end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the global discrete system
    femm = FEMMLinBeam(IntegDomain(fes, GaussRule(1, 1), B * H), material)

    K = FEMMLinBeamModule.stiffness(femm, geom0, u0, Rfield0, dchi);
    # @show K_ff = Matrix(matrix_blocked(K, nfreedofs(dchi))[:ff])
    linked = selectnode(fens; box = [h h 0 0 -d 0], inflate = h/1000)
    rb = FEMMRigidLink(
        IntegDomain(FESetP1(reshape(linked, length(linked), 1)), PointRule()), 
        linked[1], 
        Gamma)

    Krb = FEMMRigidLinkModule.stiffness(rb, geom0, u0, Rfield0, dchi)

    loaded = selectnode(fens; box = [h h 0 0 -d -d], inflate = h/1000)
    loadbdry = FESetP1(reshape(loaded, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    q = fill(0.0, 6); q[3] = -Q
    fi = ForceIntensity(q)
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve the static problem
    solve_blocked!(dchi, K + Krb, F)

    @show dchi.values[linked, :]

    if visualize
        scaling = 1e3
        dchi.values .*= scaling
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[-h -h -h];
         [h h h]]),
                    plot_nodes(fens),
                    plot_solid(fens, fes;
                               x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
                    dims = 1)
        pl = render(plots)
    end

    # # Using the stiffness coefficients of a beam in the basic configuration, we
    # # can calculate the deflection from the relationship between the shear force
    # # in the beams (there are two!) and the applied force:

    # w2 = F / (12 * E * I / h**3) / 2
    # print("Analytical deflection in the direction of the force: ", w2)
    # print("Expected tip deflection: ", [0.0, -1.25250329e-01, 0])
    # for j in m["joints"].values():
    #     print(j["displacements"])

    # plots.setup(m)
    # plots.plot_members(m)
    # plots.plot_deformations(m)
    # plots.show(m)

    # # Again using the basic stiffness, we deduce the correct magnitude of the bending moment:
    # M = 6 * h * E * I / h**3 * w2
    # print("Analytical bending moment magnitude: ", M)
    # plots.setup(m)
    # plots.plot_members(m)
    # plots.plot_bending_moments(m)
    # plots.show(m)
end # function test

function allrun()
    println("#####################################################")
    println("# test ")
    test()
end


@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

nothing
end # module
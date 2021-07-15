# disp('Simply supported square plate with center load');
module simply_supp_square_plate_conc

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellDSG3Module: FESetShellDSG3, local_frame!
using FinEtoolsFlexStructures.FEMMShellDSG3Module: FEMMShellDSG3, stiffness
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!


function doit()
    E = 30e6;
    nu = 0.3;
    force = 40;
    thickness = 0.1;
# analytical solution for the vertical deflection under the load
    analyt_sol=-0.0168;

# Mesh
    L= 10;
    n=2;
    tolerance = L/n/1000
    fens, fes = T3block(L/2,L/2,n,n);
    fens.xyz = xyz3(fens)

    sfes = FESetShellDSG3()
    accepttodelegate(fes, sfes)

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    femm = FEMMShellDSG3(IntegDomain(fes, TriRule(1), thickness), mater)

# Construct the requisite fields, geometry and displacement
# Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

# Apply EBC's
# plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 0 L/2 -Inf Inf], tolerance = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
# plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[0 L/2 0 0 -Inf Inf], tolerance = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
# simple support
    l1 = selectnode(fens; box = Float64[L/2 L/2 0 L/2 -Inf Inf], tolerance = tolerance)
    for i in [3,4,6]
        setebc!(dchi, l1, true, i)
    end
    l1 = selectnode(fens; box = Float64[0 L/2 L/2 L/2 -Inf Inf], tolerance = tolerance)
    for i in [3,5,6]
        setebc!(dchi, l1, true, i)
    end
# in-plane
    l1 = selectnode(fens; box = Float64[0 L/2 0 L/2 -Inf Inf], tolerance = tolerance)
    for i in [1, 2]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

# Assemble the system matrix
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

# # Load
# #   weight of the shell is 90 lb per square ft.
# nl=fenode_select (fens,struct ('box',[0 0 0 0  -1000  1000],'inflate',0.1));
# F = start (sysvec, get(ur, 'neqns'));
# F = assemble (F, loads(nodal_load(struct ('id',nl,'dir',[3],'magn',[-force/4])), ur));
# F = finish(F);
# # Solve
# ur = scatter_sysvec(ur, get(K,'mat')\get(F,'vec'));
# eqn=gather (ur,nl,'eqnums');
# uv=gather_sysvec(ur);
# disp(['Vertical deflection under the load: ' num2str(uv(eqn(3))) '  --  ' num2str(uv(eqn(3))/analyt_sol*100) '%'])

# # Plot
# gv=graphic_viewer;
# gv=reset (gv,[]);
# scale=100;
# u=slice(ur,(1:3),'u');
# draw(feb,gv, struct ('x', x, 'u', 0*u, 'facecolor','none'));
# # draw(feb,gv, struct ('x', x, 'u', +scale*u,'facecolor','red'));

# U=get(u,'values');
# dcm=data_colormap(struct ('range',[min(U(:,3)),max(U(:,3))], 'colormap',jet));
# colorfield=field(struct ('name', ['colorfield'], 'data',map_data(dcm, U(:,3))));
# draw(feb, gv, struct ('x', x,'u', scale*u, 'colorfield',colorfield, 'shrink',0.9));
# %
# for i=1:length(fens)
#   draw(fens(i),gv, struct ('x', x, 'u', 0*u, 'facecolor','blue'));
# end

# draw_axes(gv, struct([]));

end

end # simply_supp_square_plate_conc

simply_supp_square_plate_conc.doit()
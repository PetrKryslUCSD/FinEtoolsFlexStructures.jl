# disp('Simply supported square plate with center load');
module simply_supp_square_plate_conc

using FinEtools
using FinEtoolsFlexStructures.FESetShellDSG3Module: FESetShellDSG3
using FinEtoolsFlexStructures.FEMMShellDSG3Module: FEMMShellDSG3
using FinEtoolsDeforLinear

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
fens, fes = T3block(L/2,L/2,n,n);
fens.xyz = xyz3(fens)
fes = FESetShellDSG3(connasarray(fes); thickness = thickness)

mater = MatDeforElastIso(DeforModelRed3D, E, nu)
femm = FEMMShellDSG3(IntegDomain(fes, TriRule(1)), mater)
# # Material
# mater = mater_defor_ss_linelastic_biax (struct('E',E,'nu',nu,'type', ['stress']));
# # Finite element block
# feb = feblock_defor_ss_shell (struct ('mater',mater, 'gcells',gcells, ...
#     'integration_order',integration_order));
# # Geometry
# x = field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
# # Define the displacement/rotation field
# ur  = field(struct ('name',['ur'], 'dim', 6, 'data', zeros(get(x,'nfens'),6)));

# # Apply EBC's
# # plane of symmetry perpendicular to X
# nl=fenode_select (fens,struct ('box',[0 0   0 L/2   -10000 10000],'inflate',0.001));
# ur = set_ebc(ur, nl, nl*0+1, nl*0+1, nl*0);
# ur = set_ebc(ur, nl, nl*0+1, nl*0+5, nl*0);
# ur = set_ebc(ur, nl, nl*0+1, nl*0+6, nl*0);
# ur = apply_ebc (ur);
# # plane of symmetry perpendicular to Y
# nl=fenode_select (fens,struct ('box',[0 L/2  0 0    -10000 10000],'inflate',0.001));
# ur = set_ebc(ur, nl, nl*0+1, nl*0+2, nl*0);
# ur = set_ebc(ur, nl, nl*0+1, nl*0+4, nl*0);
# ur = set_ebc(ur, nl, nl*0+1, nl*0+6, nl*0);
# ur = apply_ebc (ur);
# # simple support
# nl=fenode_select (fens,struct ('box',[L/2 L/2  0 L/2    -10000 10000],'inflate',0.001));
# ur = set_ebc(ur, nl, nl*0+1, nl*0+3, nl*0);
# ur = set_ebc(ur, nl, nl*0+1, nl*0+4, nl*0);
# ur = set_ebc(ur, nl, nl*0+1, nl*0+6, nl*0);
# ur = apply_ebc (ur);
# nl=fenode_select (fens,struct ('box',[0 L/2   L/2 L/2   -10000 10000],'inflate',0.001));
# ur = set_ebc(ur, nl, nl*0+1, nl*0+3, nl*0);
# ur = set_ebc(ur, nl, nl*0+1, nl*0+5, nl*0);
# ur = set_ebc(ur, nl, nl*0+1, nl*0+6, nl*0);
# ur = apply_ebc (ur);
# # in-plane
# nl=fenode_select (fens,struct ('box',[0 L/2    0 L/2   -10000 10000],'inflate',0.001));
# ur = set_ebc(ur, nl, nl*0+1, nl*0+1, nl*0);
# ur = set_ebc(ur, nl, nl*0+1, nl*0+2, nl*0);
# ur = apply_ebc (ur);
# # OB(ur)
# # Number equations
# ur   = numbereqns (ur);
# # Assemble the system matrix
# K = start (sparse_sysmat, get(ur, 'neqns'));
# K = assemble (K, stiffness(feb, x, ur));
# K = finish (K);
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
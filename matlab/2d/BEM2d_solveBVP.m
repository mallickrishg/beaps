% Solve linear elastic BVP using BEM in 2-d plane strain
% boundary conditions are mixed (displacement/tractions)
% Rishav Mallick, 2023, Caltech Seismolab

clear
import geometry.*

% Elaticity parameters
G = 1;% Shear Modulus in MPa
nu = 0.25;% Poisson ratio

eM = geometry.LDhs(G,nu);
% geometry parameters(box)
Te = 1;
Ldomain = Te*2;
nbox = 40;
% parameters of ellipse
z0 = 0.6*Te;
rz = 0.3*Te;
rx = rz;
nellipse = 100;

rcv = construct_geometry(eM,Te,Ldomain,nbox,z0,rx,rz,nellipse);

figure(11),clf
plotpatch2d(rcv,rcv.xc(:,2)), hold on
quiver(rcv.xc(:,1),rcv.xc(:,2),rcv.nv(:,1),rcv.nv(:,2),0.2,'r')
axis tight equal
grid on, box on

%% apply boundary conditions
% boundary conditions are applied in the following order - left,bottom,right,internal
% 0 - Dirichlet, 1 - Neuman (set to 0)
bcindex = [0,1,0,1];

% specify BC values at each external boundary
bc_left = [-1,0];
bc_right = [1,0];

[slip_d,slip_n] = compute_bemsolution(rcv,bcindex,bc_left,bc_right);

%% compute displacements on surface
% displacement kernels and surface deformation
deltax = 0.01;
ox = linspace(-Ldomain+deltax,Ldomain-deltax,1000)';
obs = [ox,0.*ox];

[ux,uz] = compute_bemdisplacements(rcv,obs,slip_d,slip_n);

%% reference simulation
% create reference geometry object and compute solution
rcv_ref = construct_boxgeometry(eM,Te,Ldomain,nbox);
[sol_d_ref,sol_n_ref] = compute_bemsolution(rcv_ref,bcindex,bc_left,bc_right);
[ux_ref,uz_ref] = compute_bemdisplacements(rcv_ref,obs,sol_d_ref,sol_n_ref);

%% plot results
figure(10),clf
plot(ox,ux,'-','Linewidth',2), hold on
plot(ox,uz,'-','Linewidth',2)
plot(ox,ux_ref,'b--','Linewidth',2)
plot(ox,uz_ref,'r--','Linewidth',2)
axis tight, grid on, box on
xlabel('x'), ylabel('Displacement')
legend('u_x','u_z','u_x(ref)','u_z(ref)','location','northwest','box','off')
ylim([-1 1])
% xlim([-1 1]*2)
title(['r = ' num2str(rx) '; z_0 = ' num2str(z0)],'FontWeight','normal')
set(gca,'Fontsize',15,'Linewidth',1)
print(['Figures/example_r_z_uperturbation_' num2str(100*rx) '_' num2str(100*z0)],'-djpeg','-r300')

rmse = sqrt(sum((ux-ux_ref).^2 + (uz-uz_ref).^2))./length(ox);

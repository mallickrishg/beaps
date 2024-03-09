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
nbox = 10;
% parameters of ellipse
z0 = 0.5*Te;
rz = 0.3*Te;
rx = 0.3*Te;
nellipse = 100;

rcv = construct_geometry(eM,Te,Ldomain,nbox,z0,rx,rz,nellipse);

figure(11),clf
plotpatch2d(rcv,rcv.xc(:,2)), hold on
quiver(rcv.xc(:,1),rcv.xc(:,2),rcv.nv(:,1),rcv.nv(:,2),0.2,'r')
axis tight equal
grid on, box on

%% apply boundary conditions
% boundary conditions are applied in the following order - left,bottom,right,internal
% 0 - Dirichlet, 1 - Neuman, 2 - friction
bcindex = [0,1,0,2];
f_coeff = 0.6;

% specify BC values at each external boundary
bc_left = [-1,0];
bc_right = [1,0];

[slip_d,slip_n] = testing_compute_frictionbem(rcv,bcindex,bc_left,bc_right,f_coeff);

figure(1),clf
subplot(2,1,1)
plotpatch2d(rcv,slip_d)
colorbar
clim([-1 1]*max(abs(slip_d)))
axis tight equal, box on
set(gca,'Linewidth',1.5,'Fontsize',15)
subplot(2,1,2)
plotpatch2d(rcv,slip_n)
colorbar
clim([-1 1]*max(abs(slip_n)))
axis tight equal, box on
colormap("turbo")
set(gca,'Linewidth',1.5,'Fontsize',15)

%% compute tractions
[Kdd,Kdn,Knd,Knn] = geometry.computeTractionKernels(rcv,rcv);

tau_d = Kdd*slip_d + Knd*slip_n;
tau_n = Kdn*slip_d + Knn*slip_n;

figure(3),clf
subplot(2,1,1)
plotpatch2d(rcv,tau_d)
colorbar
clim([-1 1]*max(abs(tau_n)))
axis tight equal, box on
set(gca,'Linewidth',1.5,'Fontsize',15)
subplot(2,1,2)
plotpatch2d(rcv,tau_n)
colorbar
clim([-1 1]*max(abs(tau_n)))
axis tight equal, box on
colormap("bluewhitered")
set(gca,'Linewidth',1.5,'Fontsize',15)

figure(2),clf
subplot(2,1,1)
az = atan2(rcv.xc(:,2)+z0,rcv.xc(:,1));
polarplot(az(rcv.Vpl==4),slip_d(rcv.Vpl==4),'o-','LineWidth',1), hold on
polarplot(az(rcv.Vpl==4),slip_n(rcv.Vpl==4),'o-','LineWidth',1)
axis tight, box on
set(gca,'Linewidth',1.5,'Fontsize',15)
subplot(2,1,2)
polarplot(az(rcv.Vpl==4),abs(tau_d(rcv.Vpl==4)),'o-','LineWidth',1), hold on
polarplot(az(rcv.Vpl==4),abs(tau_n(rcv.Vpl==4)),'o-','LineWidth',1)
axis tight, box on
set(gca,'Linewidth',1.5,'Fontsize',15)

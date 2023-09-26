% create a bounding box and then construct a spherical cryovolcanic chamber
% use BEM to compute the perturbation in the medium due to the presence of 
% heterogeneous materials - 2d solution
% Rishav Mallick, 2023 , Caltech Seismolab

clear
import geometry.*

% elasticity material properties
Gshear = 1;% MPa
nu = 0.25;
eM = geometry.LDhs(Gshear,nu);
% geometry parameters(box)
Te = 1.0;
Ldomain = Te*2;
nbox = 40;

% parameters of ellipse
rvec = [0.1:0.05:0.4].*Te;
zvec = [0.1:0.05:0.9]*Te;
[rgrid,zgrid] = meshgrid(rvec,zvec);
index = zgrid-rgrid>=0.09 & zgrid + rgrid<=Te-0.09;
nellipse = 100;

% boundary conditions
bcindex = [0,1,0,1];
% specify BC values at each external boundary
bc_left = [-1,0];
bc_right = [1,0];

% surface displacement observations
deltax = 0.01;
ox = linspace(-Ldomain+deltax,Ldomain-deltax,1000)';
obs = [ox,0.*ox];

% create reference geometry object and compute solution
rcv_ref = construct_boxgeometry(eM,Te,Ldomain,nbox);
[sol_d_ref,sol_n_ref] = compute_bemsolution(rcv_ref,bcindex,bc_left,bc_right);
[ux_ref,uz_ref] = compute_bemdisplacements(rcv_ref,obs,sol_d_ref,sol_n_ref);

%% loop through various r,z0 values for cryovolcano chamber 
UXvec = nan(length(ux_ref(:)),length(rgrid(:)));
UZvec = nan(length(ux_ref(:)),length(rgrid(:)));

for iter = 1:length(rgrid(:))
    if index(iter)==1
        r = rgrid(iter);% dz of volcano chamber
        z0 = zgrid(iter);% depth of chamber        
        % construct bounding box and volcano meshes
        
        rcv = construct_geometry(eM,Te,Ldomain,nbox,z0,r,r,nellipse);                              
        % compute solution and surface displacements
        [sol_d,sol_n] = compute_bemsolution(rcv,bcindex,bc_left,bc_right);                
        
        tic
        [ux,uz] = compute_bemdisplacements(rcv,obs,sol_d,sol_n);        
        toc
        
        UXvec(:,iter) = ux;
        UZvec(:,iter) = uz;
        
        disp(['Simulation number:' num2str(iter)])
    end
end

%% plot specific results
figure(2),clf
subplot(211)
plot(ox,UXvec), hold on
plot(ox,ux_ref,'r-','Linewidth',2)
axis tight
ylim([-1 1])

subplot(212)
plot(ox,UZvec), hold on
plot(ox,uz_ref,'r-','Linewidth',2)
axis tight
ylim([-1 1])

figure(10),clf
imagesc(UZvec)
caxis([-1 1]*0.2)

figure(11),clf
toplot = sqrt(sum((UZvec-uz_ref).^2 + (UXvec-ux_ref).^2))./length(ox);
scatter(rgrid(:),zgrid(:),200,toplot,'s','filled','LineWidth',1,'MarkerEdgeColor','k')
xlabel('r'),ylabel('z_0')
xlim([0 0.5])
ylim([0 1])
box on
caxis([0 1]*5e-3)
cb=colorbar;cb.Label.String = '\Deltau';cb.LineWidth=1;
set(gca,'Fontsize',20,'Linewidth',1.5)
print('Figures/vary_r_z_uperturbation','-djpeg','-r300')
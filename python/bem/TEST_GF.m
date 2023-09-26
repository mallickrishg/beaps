%  testing greens functions
% testing how displacement greens functions work
% Rishav Mallick, Caltech, 2022

clear
addpath ~/Dropbox/scripts/utils/

% shear modulus
G = 30e3;
nu = 0.25;

% grid points
nx2 = 100;
nx3 = nx2;

x2 = linspace(-50e3,50e3,nx2);
x3 = linspace(-50e3,50e3,nx3);
[X2,X3] = meshgrid(x2,x3);

% source properties
Y2 = 0e3;
Y3 = 0e3;
W = 20e3;% half-length of crack
dip = 45;
slip = 1;
open = 0;

%% compute displacements in full-space
ng = 13;
[Disp] = LDdispFS(X2,X3,Y2,Y3,W/2,-deg2rad(dip),slip,open,nu);
ue = Disp(:,1);uz = Disp(:,2);


figure(11),clf
subplot(211)
imagesc(x2./1e3,x3./1e3,reshape(ue,nx3,nx2)), hold on
contour(x2./1e3,x3./1e3,reshape(ue,nx3,nx2),[-1:.1:1].*0.25,'k-')
quiver(X2(1:ng:end)./1e3,X3(1:ng:end)./1e3,ue(1:ng:end)',uz(1:ng:end)','k-','Linewidth',1)
plot(x2./1e3,0.*x2,'k-','LineWidth',1)
axis tight equal, box on
cb = colorbar;cb.Label.String='u_x';
caxis([-1,1].*0.5)
set(gca,'Fontsize',15,'YDir','normal','LineWidth',2)

subplot(212)
imagesc(x2./1e3,x3./1e3,reshape(uz,nx3,nx2)), hold on
contour(x2./1e3,x3./1e3,reshape(uz,nx3,nx2),[-1:.1:1].*0.25,'k-')
plot(x2./1e3,0.*x2,'k-','LineWidth',1)
axis tight equal, box on
cb = colorbar;cb.Label.String='u_z';
caxis([-1,1].*0.5)
colormap(bluewhitered(40))
set(gca,'Fontsize',15,'YDir','normal','LineWidth',2)

%% compute Stress
[Stress] = LDstressFS(X2,X3,Y2,Y3,W/2,-deg2rad(dip),slip,open,nu,2*G*(1+nu));
sxx = Stress(:,1);syy = Stress(:,2);sxy = Stress(:,3);

figure(12),clf
imagesc(x2./1e3,x3./1e3,reshape(sxy,nx3,nx2)),
axis tight equal
caxis([-1 1])


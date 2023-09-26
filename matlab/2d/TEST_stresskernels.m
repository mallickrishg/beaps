
clear

x = linspace(-10e3,10e3,100);

% shear modulus
G = 30e3;
nu = 0.25;

% grid points
nx2 = 11;

x2 = linspace(-10e3,10e3,nx2);
x3 = 0.*x2;

% source properties
Y2 = 0e3;
Y3 = 0e3;
W = 2e3;% length of crack
dip = 0;
slip = 1;
open = 0;

[Stress] = LDstressFS(x2,x3,Y2,Y3,W/2,-deg2rad(dip),slip,open,nu,2*G*(1+nu));
sxx = Stress(:,1);szz = Stress(:,2);sxz = Stress(:,3);

% plot figures
figure(1),clf
subplot(311)
plot(x2./1e3,sxx,'.-','Linewidth',2)
subplot(312)
plot(x2./1e3,szz,'.-','Linewidth',2)
subplot(313)
plot(x2./1e3,sxz,'o-','Linewidth',2)
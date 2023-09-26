function rcv = construct_geometry(eM,Te,Le,nbox,z0,rx,rz,nellipse)

% create a box of dimensions 2Le x Te
patchfname = 'box2d.seg';

fileID = fopen(patchfname,'w');
fprintf(fileID,'%s\n',...
    '# n  Vpl  X  Z(+)  Width  Dip  W0  qW');
% left face
fprintf(fileID,'%d %.2f %.12f %.12f %.12f %.12f %.12f %d\n',...
    1, 1, -Le, 0, Te, 90, Te/nbox, 1);
% bottom face
fprintf(fileID,'%d %.2f %.12f %.12f %.12f %.12f %.12f %d\n',...
    1, 2, -Le, Te, 2*Le, 0, Te/nbox, 1);
% right face
fprintf(fileID,'%d %.2f %.12f %.12f %.12f %.12f %.12f %d\n',...
    1, 3, Le, Te, Te, 270, Te/nbox, 1);

% mesh an ellpise at (0,z0) with 2 radii specified
theta = -linspace(0,-2*pi,nellipse);
x = rx*cos(theta);
z = rz*sin(theta) + z0;

for i = 1:length(x)-1
    w = sqrt((x(i+1)-x(i))^2 + (z(i+1)-z(i))^2);
    dip = atan2d(z(i+1)-z(i),x(i+1)-x(i));
    fprintf(fileID,'%d %.2f %.12f %.12f %.12f %.12f %.12f %d\n',...
        1, 4, x(i), z(i), w, dip, w, 1);
end

fclose(fileID);
rcv = geometry.receiver(patchfname,eM);

end
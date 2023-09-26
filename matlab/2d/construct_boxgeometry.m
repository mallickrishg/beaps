function rcv = construct_boxgeometry(eM,Te,Le,nbox)

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

fclose(fileID);
rcv = geometry.receiver(patchfname,eM);

end
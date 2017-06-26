clear;
close all;
%% 
% nz = 222;
% nx = 222;
% ny = 41;
% dx = 1;
% dz = 1;

% xmin = -49; xmax = xmin + (nx-1) * dx;
% zmin = -5000; zmax = zmin + (nz-1) * dz;
% ForwardX = xmin:dx:xmax;
% ForwardZ = zmin:dz:zmax;
% model = load('./Par/model3D.txt'); 
% model = reshape(model(:,3), nz, nx, ny);
% forward grids
dx = 1; xmin = 0; xmax = 221;
dy = 1; ymin = 0; ymax = 221;
dz = 1; zmin = -36.0; zmax = 4;
ForwardX = xmin:dx:xmax;
ForwardY = ymin:dy:ymax;
ForwardZ = zmin:dz:zmax;
%%
!../main
%%
rays = load('../Out/raypoints.txt');
% imagesc(ForwardX, ForwardZ, model);
% colormap('jet');5
% hold on;
figure;
scatter3(rays(:,3), rays(:,4), rays(:,5), 'r', '.');
% plot3(rays(:,3),rays(:,4),rays(:,5),'r.-');
box on; axis image; xlabel('x ');ylabel('y ');zlabel('z ');
% axis([ForwardX(1) ForwardX(end) ForwardY(1) ForwardY(end) ForwardZ(1) ForwardZ(end)]);
hold on;
ax = [rays(1,3), rays(end,3)]; ay = [rays(1,4), rays(end,4)]; az = [rays(1,5), rays(end,5)];
plot3(ax, ay, az);
view(-112,10);
gcf.PaperPositionMode = 'auto';
print('rays','-djpeg','-r300');
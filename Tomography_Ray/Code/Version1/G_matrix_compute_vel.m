function G = G_matrix_compute_vel(Vp, raypoints, nz, oz, dz, nx, ox, dx)
%%

% hold on;
% scatter(mid_raypoints(:,2), mid_raypoints(:,1),'.');
% xlabel('x');
% ylabel('z');
% axis ij;

G = zeros(nz,nx);

for i=1:size(raypoints,1)-1
    mid_raypoints(i,:) = (raypoints(i,:) + raypoints(i+1,:)) /2;
end

for i=20:size(mid_raypoints,1)-20
    zIndex_up = floor((mid_raypoints(i,1) - oz)/dz+1);
    zIndex_down = zIndex_up + 1;
    xIndex_left = floor((mid_raypoints(i,2) - ox)/dx+1);
    xIndex_right = xIndex_left + 1;
%     wz = zIndex_down - mid_raypoints(i,1);
%     wx = xIndex_right - mid_raypoints(i,2);
    wz = (oz + (zIndex_down-1)*dz - mid_raypoints(i,1))/dz;
    wx = (ox + (xIndex_right-1)*dx - mid_raypoints(i,2))/dx;
%     if (i==size(mid_raypoints,1))
        dh = sqrt((raypoints(i+1,1)-raypoints(i,1))^2 + (raypoints(i+1,2)-raypoints(i,2))^2);
%     end

    vel = Vp(zIndex_up, xIndex_left)*wz*wx ...
        + Vp(zIndex_up, xIndex_right)*wz*(1-wx) ...
        + Vp(zIndex_down, xIndex_left)*(1-wz)*wx ...
        + Vp(zIndex_down, xIndex_right)*(1-wz)*(1-wx);
    
    G(zIndex_up, xIndex_left) = G(zIndex_up, xIndex_left) - wz*wx*dh/vel^2;
    G(zIndex_up, xIndex_right) = G(zIndex_up, xIndex_right) - wz*(1-wx)*dh/vel^2;
    G(zIndex_down, xIndex_left) = G(zIndex_down, xIndex_left) - (1-wz)*wx*dh/vel^2;
    G(zIndex_down, xIndex_right) = G(zIndex_down, xIndex_right) - (1-wz)*(1-wx)*dh/vel^2;
    
    
end

% imagesc(G);


G = G(:);
G = G';
end

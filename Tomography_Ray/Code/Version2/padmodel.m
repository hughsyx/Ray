function [ velPadded ] = padmodel( velOriginal, PADEDGE )
%pad edges of velocity model for computing
%   PADEDGE is the edge grid number for padding
%   Assume: velOriginal(nz,nx,ny)
[nz,nx,ny] = size(velOriginal);
velAve = mean(velOriginal(:));
if (ny==1) %2D
    velPadded = velAve*ones(nz+2*PADEDGE,nx+2*PADEDGE);
    % pad x
    velPadded(PADEDGE+1:PADEDGE+nz,:) = ...
    [repmat(velOriginal(:,1),1,PADEDGE) velOriginal repmat(velOriginal(:,end),1,PADEDGE)];
    % pad z
    velPadded(1:PADEDGE,:) = repmat(velPadded(PADEDGE+1,:),PADEDGE,1);
    velPadded(PADEDGE+nz+1:2*PADEDGE+nz,:) = repmat(velPadded(PADEDGE+nz,:),PADEDGE,1);
else %3D
    velPadded = velAve*ones(nz+2*PADEDGE,nx+2*PADEDGE,ny+2*PADEDGE);
    % pad x
    velPadded(PADEDGE+1:PADEDGE+nz,:,PADEDGE+1:PADEDGE+ny) = ...
    [repmat(velOriginal(:,1,:),1,PADEDGE,1) velOriginal repmat(velOriginal(:,end,:),1,PADEDGE,1)];
    % pad z
    velPadded(1:PADEDGE,:,PADEDGE+1:PADEDGE+ny) = ...
        repmat(velPadded(PADEDGE+1,:,PADEDGE+1:PADEDGE+ny),PADEDGE,1,1);
    velPadded(PADEDGE+nz+1:2*PADEDGE+nz,:,PADEDGE+1:PADEDGE+ny) = ...
        repmat(velPadded(PADEDGE+nz,:,PADEDGE+1:PADEDGE+ny),PADEDGE,1,1);
    % pad y
    velPadded(:,:,1:PADEDGE) = repmat(velPadded(:,:,PADEDGE+1),1,1,PADEDGE);
    velPadded(:,:,PADEDGE+ny+1:2*PADEDGE+ny) = ...
        repmat(velPadded(:,:,PADEDGE+ny),1,1,PADEDGE);
end

end


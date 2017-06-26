%%%%  Crosswell 2D Ray Tracing Program   %%%%
%%%% by Shaoyu Lu, Nov, 2013 at Stanford %%%%
function [DAOtime,G,timeTableAllShots] = DAOtimeG(velPadded,padForwardX,padForwardY,padForwardZ,...
    PADEDGE,InvertX,InvertY,InvertZ,sou_location,rec_location,dxz)
%% forward grids
xmin = min(padForwardX);
dx = abs(padForwardX(2)-padForwardX(1));
zmin = min(padForwardZ);
dz = abs(padForwardZ(2)-padForwardZ(1));
%% inverse grids
invnx = length(InvertX);
invnz = length(InvertZ);
%% source and receiver configuration
zsrcAll = sou_location(:,end);
nsrc = length(zsrcAll);
xsrcAll = sou_location(:,1);
zrecAll = rec_location(:,end);
nrec = length(zrecAll);
xrecAll = rec_location(:,1);
%% compute traveltime map for each shot
DAOtime=zeros(nrec,nsrc);
timeTableAllShots = cell(nsrc,1);
for ishot=1:nsrc
    xsrc = xsrcAll(ishot); 
    ixsrc = round((xsrc-xmin)/dx)+1;
    zsrc = zsrcAll(ishot); 
    izsrc = round((zsrc-zmin)/dz)+1;
    isrcLocZX = [izsrc,ixsrc]'; %padded velocity causes index shift
    fprintf(['\nComputing traveltime table for shot#',num2str(ishot),' ...\n']);
    tic
    timeTableAllShots{ishot} = msfm(velPadded, isrcLocZX, true, false)*dxz;
    toc
    for irec = 1:nrec
        xrec = xrecAll(irec); 
        ixrec = round((xrec-xmin)/dx)+1;
        zrec = zrecAll(irec); 
        izrec = round((zrec-zmin)/dz)+1;
        DAOtime(irec,ishot) = timeTableAllShots{ishot}(izrec,ixrec);
    end
end
%% compute inversion matrix G by ray-tracing based on traveltime maps
tic
G = zeros(nsrc*nrec,(invnx-1)*(invnz-1));
% shot gather rays
ysrc=0; yrec=0; %2D
for ishot=1:nsrc
    xsrc = xsrcAll(ishot); 
    zsrc = zsrcAll(ishot);
    for irec=1:nrec
        xrec = xrecAll(irec);
        zrec = zrecAll(irec);
        iray = (ishot-1)*nrec+irec;
        %fprintf(['\nComputing raytracing for shot#',num2str(ishot),' and rec#',num2str(irec),' ...\n']);
        %tic
        [G(iray,:),raypoints1] = raytracing(timeTableAllShots{ishot},padForwardX,padForwardY,...
            padForwardZ,InvertX,InvertY,InvertZ,xsrc,ysrc,zsrc,xrec,yrec,zrec);
        %toc
        %% plot rays
        %{
            G1 = G(iray,:);
            figure(100); subplot(121);
            imagesc(padForwardX(1+PADEDGE:end-PADEDGE),padForwardZ(1+PADEDGE:end-PADEDGE),...
                timeTableAllShots{ishot}(1+PADEDGE:end-PADEDGE,1+PADEDGE:end-PADEDGE));
            H=colorbar;set(get(H,'Title'),'string','s ');axis image;
            hold on; plot(raypoints1(1,:),raypoints1(3,:),'r.-');
            midInvertX = (InvertX(1:end-1)+InvertX(2:end))/2;% middle of cells
            midInvertZ = (InvertZ(1:end-1)+InvertZ(2:end))/2;% middle of cells
            subplot(122); imagesc(midInvertX,midInvertZ,reshape(G1,invnz-1,invnx-1)); 
            axis image; colorbar;
            pause(0.2);
        %}
    end
end
TimeG = toc;
fprintf(['\nTotal time for computing G matrix is ',num2str(TimeG),' sec\n']);
%{
figure;
subplot(131);
imagesc(padForwardX(1+PADEDGE:end-PADEDGE),padForwardZ(1+PADEDGE:end-PADEDGE),...
    velPadded(1+PADEDGE:end-PADEDGE,1+PADEDGE:end-PADEDGE));
xlabel('X (m) ');ylabel('Z (m) ');title('Current Velocity Model ');
H=colorbar;set(get(H,'Title'),'string','(m/s) ');axis image;
subplot(132);
imagesc(G);axis image;colorbar;
xlabel('pixel number ');ylabel('ray number');title('Inversion Matrix G ');
subplot(133);
imagesc(sou_location(:,2),rec_location(:,2),DAOtime);H=colorbar;
set(get(H,'Title'),'string','(s) ');axis image;
xlabel('Source depth (m) ');ylabel('Receiver depth (m) ');
title('Direct Arrival Traveltime Computed by Eikonal Solver ','FontWeight','Bold');
%}
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Direct-Arrival Tomographic Inversion                    %
%                        Matlab calling C++ code                          %
%                       Shaoyu Lu, SWP, Nov.2013                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
tstart = cputime;
Tomo_path = '/Users/Hugh/Documents/RAY/Tomography_Ray';

%% Pre-tomo preparation: 

% prepare:

%    Generate_vel_1circle.m
%    --This mfile:
%       1) generate the true velocity model and source/receiver geometry
%       2) generate the synthetic direct-arrivals using Eikonal solver

load([Tomo_path, '/Models/marmousi_vp.mat']);
load([Tomo_path, '/Models/marmousi_vp_smoothed.mat']);
Vtrue = marmousi_vp;
Vsmoothed = marmousi_vp_smoothed;
clear('marmousi_vp', 'marmousi_vp_smoothed');
Vtrue = Vtrue(21:120, 101:200);
Vtrue = modelSmooth(Vtrue, 3);
Vsmoothed = Vsmoothed(21:120, 101:200);
Strue = 1./Vtrue;

%% Generate starting model, inversion grids, and traveltime0, G matrix
% read source and receiver information
vecZSource = 5:2:95;
vecXSource = 5*ones(1,length(vecZSource));
vecZReceiver = 5:2:95;
vecXReceiver = 95*ones(1,length(vecZReceiver));
sCnt = length(vecZSource);
rCnt = length(vecZReceiver);
sou_location = [vecXSource' vecZSource'];
rec_location = [vecXReceiver' vecZReceiver'];

% 1. Define inversion grids and pad forward model
[nz, nx] = size(Vtrue);
ny = 1;
deltax = 1;
deltay = 1;
deltaz = 1;
xmin = 1;
ymin = 1;
zmin = 1;
ForwardX = xmin:deltax:xmin+deltax*(nx-1);
ForwardZ = zmin:deltaz:zmin+deltaz*(nz-1);
% inverse grids
%%%%%% change this %%%%%%%%%%
invnz = 21; invnx = 21;
% invnz = sCnt;invnx = round(nx/nz*invnz)+1;
% invnz = 100;invnx = round(nx/nz*invnz)+1;
invny = 1;
k0 = 2;% the iteration time for non-regularized LSQR inversion
k1 = 2;% the iteration time for regularized LSQR inversion
smoothZtoX = 1;% smooth strength along z to along x
Maxiter = 20; % how many iteration do you want to make?
stopRRMS = 0.001; % for how much relative change to stop
plotLcurve = 'no'; % or yes to plot the L curve
plotEveryInversion = 'no'; % or yes to plot every inversion result
saveResults = 'no'; % save the results?
constVel0 = 'yes';% constant velocity starting model?
bound_dS = 'yes'; bound = 0.05; % bound update? the default is 5%
smooth_dS = 'yes';% smooth update?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin = min(ForwardX); xmax = max(ForwardX);
zmin = min(ForwardZ); zmax = max(ForwardZ);
dx = abs(ForwardX(2)-ForwardX(1));
dz = abs(ForwardZ(2)-ForwardZ(1));
midInvertX = linspace(xmin,xmax,invnx-1);% make rays inside cells
InvertY = 0;
midInvertZ = linspace(zmin,zmax,invnz-1);% make rays inside cells
invdx = midInvertX(2) - midInvertX(1);
invdy = 0;
invdz = midInvertZ(2) - midInvertZ(1);
InvertX = linspace(xmin - invdx/2,xmax + invdx/2,invnx);% middle of cells
InvertZ = linspace(zmin - invdz/2,zmax + invdz/2,invnz);% middle of cells

% define a padding edge for computing
PADEDGE = max([2,round(invdx/dx),round(invdz/dz)]);
padForwardX = xmin-PADEDGE*dx:dx:xmax+PADEDGE*dx;
padForwardY = 0;%2D
padForwardZ = zmin-PADEDGE*dz:dz:zmax+PADEDGE*dz;

% True traveltime computation
Vtrue_pad = padmodel(Vtrue, PADEDGE );
 [DAOtime_obs, ~] = DAOtimeG(Vtrue_pad,padForwardX,padForwardY,...
        padForwardZ,PADEDGE,InvertX,InvertY,InvertZ,sou_location,rec_location,dx);
% 2. generate 0-D or 1-D starting model and inversion grids
if strcmp(constVel0,'yes') % const velocity starting model
    [slow0,vel0] = build0Dvel0(sou_location,rec_location,...
    ForwardX,ForwardZ,DAOtime_obs);
else % 1D starting model
    [slow0,vel0] = build1Dvel0(sou_location,rec_location,isou_location,...
    ForwardX,ForwardZ,DAOtime_obs);
end
S{1}=slow0; V{1}=vel0; 
S{1} = 1./Vsmoothed; V{1} = Vsmoothed;
Vcolorscale = [min(Vtrue(:)) max(Vtrue(:))];
Scolorscale = [min(Strue(:)) max(Strue(:))];


%% Do inversion iteratively
for iter=1:Maxiter
    iter
    % prepare for interpolation
    [xxi,zzi]=meshgrid(midInvertX,midInvertZ);
    [xxt,zzt]=meshgrid(ForwardX,ForwardZ);
    % 3. compute the direct-arrival traveltime for current model
    % pad velocity model for computing traveltime
    Vpad = padmodel( V{iter}, PADEDGE );
    % compute direct-arrival traveltime, G matrix and timeTable for all shots
    [DAOtime_cal,G] = DAOtimeG(Vpad,padForwardX,padForwardY,...
        padForwardZ,PADEDGE,InvertX,InvertY,InvertZ,sou_location,rec_location,dx);
    % compute traveltime difference 
    % (traveltime(irec,ishot) -> one column, shot gather, consistent with G)
    dTraveltime{iter} = reshape(DAOtime_obs-DAOtime_cal,sCnt*rCnt,1);
    RMS(iter) = norm(dTraveltime{iter});
    % 4. do inversion to compute slowness update
%     [dS{iter},lambda_c(iter),InvTime(iter),RMS(iter+1)]=LSQRinversionS_D2(...
%         plotLcurve,G,dTraveltime{iter},k0,k1,midInvertX,midInvertZ,...
%         sCnt,rCnt,smoothZtoX);
        [dS{iter},lambda_c(iter),InvTime(iter),RMS(iter+1)]=LSQRinversionS(...
        G,dTraveltime{iter},k0,k1,midInvertX,midInvertZ,...
        sCnt,rCnt,smoothZtoX);
    % 5. interp the coarse grid inversion result with fine grids
    dS{iter}=reshape(dS{iter},invnz-1,invnx-1);
    dS_f{iter}=interp2(xxi,zzi,dS{iter},xxt,zzt,'spline');
    %% smooth dS update
    if strcmp(smooth_dS,'yes')
%         dS_f{iter}=smoothn(dS_f{iter},1000);
        dS_f{iter} = modelSmooth(dS_f{iter}, 3);
        dS{iter}=interp2(xxt,zzt,dS_f{iter},xxi,zzi,'spline');
        %dS{iter}=smoothn(dS{iter});
        RMS(iter+1)=norm(G*reshape(dS{iter},(invnz-1)*(invnx-1),1)-dTraveltime{iter});
    end
    dS_per=dS_f{iter}./S{iter};
    %% bound the pixel-wise dS update 
    if strcmp(bound_dS,'yes')
        index = find(abs(dS_per)>bound);
        dS_f{iter}(index)=S{iter}(index).*sign(dS_per(index))*bound;
        dS_per=dS_f{iter}./S{iter};
    end
    
    %%
    dS_per_max(iter)=max(abs(dS_per(:)))*100;
    S{iter+1}=dS_f{iter}+S{iter};
    V{iter+1}=1./S{iter+1};   
    % 6. plot the updated velocity
    if strcmp(plotEveryInversion,'yes')
    figure;
    subplot(221);set(gca,'fontsize',12);
    imagesc(ForwardX,ForwardZ,vtrue,Vcolorscale);
    H=colorbar;set(get(H,'Title'),'string','m/s '); hold on;
    plot(sou_location(:,1),sou_location(:,end),'r*','markersize',8);
    hold on;plot(rec_location(:,1),rec_location(:,end),'g<','markersize',8);
    title('True Velocity ','FontSize',18,'fontweight','b');
    xlabel('Distance (m) ','FontSize',18);
    ylabel('Depth (m) ','FontSize',18);axis image;
    subplot(222);set(gca,'fontsize',12);
    imagesc(ForwardX,ForwardZ,V{iter+1},Vcolorscale);
    H=colorbar;set(get(H,'Title'),'string','m/s '); hold on;
    plot(sou_location(:,1),sou_location(:,end),'r*','markersize',8);
    hold on;plot(rec_location(:,1),rec_location(:,end),'g<','markersize',8);
    title(['Updated Velocity by ',num2str(iter),' iteration'],...
        'FontSize',18,'fontweight','b');
    xlabel('Distance (m) ','FontSize',18);
    ylabel('Depth (m) ','FontSize',18);axis image;
    subplot(223);set(gca,'fontsize',12);
    plot([0:iter],RMS,'b+-');axis square;
    xlabel('Iteration times ','FontSize',18);
    ylabel('Traveltime RMS error ','FontSize',18);
    title({'Traveltime RMS error ','change with iterations '},'FontSize',18,'fontweight','b');
    subplot(224);set(gca,'fontsize',12);
    plot([1:iter],dS_per_max,'b+-');axis square;
    xlabel('Iteration times ','FontSize',18);
    ylabel('Absolute Percentage Change (%) ','FontSize',18);
    title({'Maximum Pixel-wise Slowness Update ','change with iterations '},'FontSize',18,'fontweight','b');
    drawnow;
    end
    % 7. see if stop the iteration or not
    if abs(RMS(iter+1)-RMS(iter))/RMS(iter) < stopRRMS 
        break;
    end
end
TotalTime = cputime - tstart;
fprintf(['\nThe total time for tomography and ploting is ',num2str(TotalTime),'s\n']);

%% save results
% if strcmp(saveResults,'yes')
% if strcmp(constVel0,'yes') % const velocity starting model
%     save ([foldername,'/tomoresults_0Dvel0_finalState.mat'],'V','S','RMS',...
%         'dS_per_max','ForwardX','ForwardZ','InvertX','InvertZ','vtrue',...
%         'Vcolorscale','TotalTime','DAOtime_cal','DAOtime_obs','G');
% else % 1D starting model
%     save ([foldername,'/tomoresults_1Dvel0_finalState.mat'],'V','S','RMS',...
%         'dS_per_max','ForwardX','ForwardZ','InvertX','InvertZ','vtrue',...
%         'Vcolorscale','TotalTime','DAOtime_cal','DAOtime_obs','G');
% end
% end
% figure; 
%     subplot(231);set(gca,'fontsize',12);
%     imagesc(ForwardX,ForwardZ,Strue,Scolorscale); colormap(jet);
%     H=colorbar;set(get(H,'Title'),'string','s/m '); hold on;
%     plot(sou_location(:,1),sou_location(:,end),'r*','markersize',8);
%     hold on;plot(rec_location(:,1),rec_location(:,end),'g<','markersize',8);
%     title('True Slowness ','FontSize',16,'fontweight','b');
%     xlabel('Distance (m) ','FontSize',14);
%     ylabel('Depth (m) ','FontSize',14);axis image;
%     subplot(232);set(gca,'fontsize',12);
%     imagesc(ForwardX,ForwardZ,1./V{1},Scolorscale); colormap(jet);
%     H=colorbar;set(get(H,'Title'),'string','s/m '); 
%     hold on;plot(sou_location(:,1),sou_location(:,end),'k*','markersize',8);
%     hold on;plot(rec_location(:,1),rec_location(:,end),'k<','markersize',8);
%     title('Starting Slowness Model ','FontSize',16,'fontweight','b');
%     xlabel('Distance (m) ','FontSize',14);
%     ylabel('Depth (m) ','FontSize',14);axis image;
%     subplot(233);set(gca,'fontsize',12);
%     imagesc(ForwardX,ForwardZ,1./V{end},Scolorscale); colormap(jet);
%     H=colorbar;set(get(H,'Title'),'string','s/m '); hold on;
%     plot(sou_location(:,1),sou_location(:,2),'r*','markersize',8);
%     hold on;plot(rec_location(:,1),rec_location(:,end),'g<','markersize',8);
%     title({'Inverted Slowness ',['by ',num2str(iter),' iteration']},...
%         'FontSize',16,'fontweight','b');
%     xlabel('Distance (m) ','FontSize',14);
%     ylabel('Depth (m) ','FontSize',14);axis image;
%     subplot(234);set(gca,'fontsize',12);
%     plot([0:iter],RMS,'b+-');axis square;
%     xlabel('Iteration times ','FontSize',14);
%     ylabel('Traveltime RMS error ','FontSize',14);
%     title({'Traveltime RMS error ','change with iterations '},'FontSize',16,'fontweight','b');
%     subplot(235);set(gca,'fontsize',12);
%     plot([1:iter],dS_per_max,'b+-');axis square;
%     xlabel('Iteration times ','FontSize',14);
%     ylabel('Absolute Percentage (%) ','FontSize',14);
%     title({'Maximum Pixel-wise Slowness Update ','Changing with iterations '},...
%         'FontSize',16,'fontweight','b');
%     subplot(236);set(gca,'fontsize',12);
%     imagesc(sou_location(:,end),rec_location(:,end),DAOtime_obs-DAOtime_cal);axis image;
%     H=colorbar;set(get(H,'Title'),'string','s ');
%     xlabel('Source Depth (m) ','FontSize',14);
%     ylabel('Receiver Depth (m) ','FontSize',14);
%     title('Traveltime Residual (Obs - Cal) ','FontSize',16,'fontweight','b');
%     drawnow;

%%
figure;
pcolor(Vtrue);
shading interp;
colormap(jet);
colorbar;
axis ij;

figure;
pcolor(V{1});
shading interp;
colormap(jet);
colorbar;
axis ij;

figure;
pcolor(V{2});
shading interp;
colormap(jet);
colorbar;
axis ij;

% figure;
% pcolor(dS);
% shading interp;
% colormap(jet);
% colorbar;
% axis ij;

figure;
plot(RMS);
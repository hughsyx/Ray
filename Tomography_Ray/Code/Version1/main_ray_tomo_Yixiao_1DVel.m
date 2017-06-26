%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              1 Dimentional Tomography
%                               Yixiao Sheng, Stanford
%                                           May 3, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
% tic;

%% set paths
Tomo_path = '/Users/Hugh/Documents/RAY/Tomography_Ray';
RT_path = '/Users/Hugh/Documents/RAY/Eikonal';

%% cr0ss-att test
load([Tomo_path, '/Models/cross_well_att_vp.mat']);
load([Tomo_path, '/Models/cross_well_att_0Dvp.mat']);
load([Tomo_path, '/Models/cross_well_att_T_obs_new.mat']);
Vtrue = cross_well_att_vp;
Vsm = V0D;
T_obs = t_esti_model2';
clear('cross_well_att_vp','V0D','T_obs_my');
% Vtrue = modelSmooth(Vtrue, 3);

%% block test
% Vp_smoothed = 3000*ones(100, 100);
% Vp = Vp_smoothed;
% Vp(40:60, 40:60) = 3400;

%% layer test
% Vsm = 3000 * ones(100, 100);
% Vtrue = Vsm;
% Vtrue(20:40, :) = 2500;
% Vtrue(60:80, :) = 3200;

%% set parameter files
nDim = 2;
bound = 0.05;
[nz, nx] = size(Vtrue);
ny = 0;
deltax = 1;
deltay = 1;
deltaz = 1;
xmin = 1;
ymin = 1;
zmin = 1;

nz_G=nz; nx_G=2; % parameter for inversion
deltax_G = 139; deltaz_G=deltaz;
% vecZSource = 100:100:5300;
% vecXSource = 100*ones(1,length(vecZSource));
% vecZReceiver = 200:100:5200;
% vecXReceiver = 600*ones(1,length(vecZReceiver));
vecZSource = 10*2:10*2:530*2;
vecXSource = 10*2*ones(1,length(vecZSource));
vecZReceiver = 20*2:10*2:520*2;
vecXReceiver = 60*2*ones(1,length(vecZReceiver));

nTer = 10;

%% create parameter and model files
writeParfile(RT_path, nDim, nx, ny, nz, deltax, deltay, deltaz, xmin, ymin, zmin, ...
    vecZSource, vecXSource, vecZReceiver, vecXReceiver);

writeModel(RT_path, modelSmooth(Vtrue, 15));

%% run ray tracing code
cd(RT_path);
! ./main
cd([Tomo_path, '/Code/Version1']);

%% compute rays, ray traveltime, and the G matrix
%[T_obs, Raypoints, ~] = readRaypoints_compGmat(RT_path, nz, ...
%    zmin, deltaz, nx, xmin, deltax);

%% plot the rays
% figure;
% pcolor(Vtrue); shading interp; colormap(jet); axis ij;
% hold on;
% for i = 1:size(Raypoints, 1)
%     plot(Raypoints{i}(:, 2), Raypoints{i}(:, 1), 'k');
%     hold on;%G_matrix
% end
% toc;

%% inversion
tic;
V_collec = zeros(nz, nx, nTer + 1);
V_collec(:, :, 1) = Vsm;%Vsm;
S = mean(1./Vsm,2);
Residual = zeros(nTer, 1);
for iTer = 1:nTer;
    disp(['Iteration ', num2str(iTer)]);
    writeModel(RT_path, modelSmooth(V_collec(:, :, iTer), 15));
    
    cd(RT_path);
    ! ./main
    cd([Tomo_path, '/Code/Version1']);
    
    [T_syn, Raypoints, G_matrix] = readRaypoints_compGmat(RT_path, nz_G, ...
        zmin, deltaz_G, nx_G, xmin, deltax_G);
    
    dT = T_obs - T_syn;
    
    Residual(iTer) = norm(dT);
    
    lambda = max(max(abs(G_matrix)));
%     dS_per = 1;
%     while(max(max(abs(dS_per))) >= bound)

        cvx_begin
        cvx_quiet(true)
        variable dS(nz_G, nx_G)
        minimize(norm(dT - G_matrix * dS(:)) + lambda * norm(dS(:)))
        cvx_end
        
%         tic;
%         disp('Total variation:');
%         cvx_begin
%         cvx_quiet(true)
%         variable dS_tv(nz, nx)
%         minimize(tv(dS_tv))
%         subject to
%         dT == G_matrix * dS_tv(:)
%         cvx_end
%         toc;

        
%         Zero_vector = zeros(size(G_matrix,2), 1);
%         I_matrix = 7 * lambda*eye(size(G_matrix,2));
%         G2 = [G_matrix; I_matrix];
%         d2 = [dT; Zero_vector];
%         B1 = lsqr(G2, d2, 1e-06, 10000);
%         dS  = reshape(B1, nz, nx);
        
        %dS = modelSmooth(dS, 3);
        dS_per = repmat(mean(dS,2)./S,1,nx);
%         max(max(abs(dS_per)))
%         lambda  = lambda * 1.2;
%         lambda
%     end
    
       index = find(abs(dS_per)>bound);
       dS_per(index) = sign(dS_per(index)) * bound;
       V_collec(:, :, iTer + 1) = V_collec(:, :, iTer)./(1 + dS_per);
       S = 1./V_collec(:, :, iTer + 1);
       S = mean(S,2);
end
toc;

%% Total variation
% tic;
% disp('Total variation:');
% cvx_begin
% cvx_quiet(true)
%     variable dS_tv(nz, nx)
%     minimize(tv(dS_tv))
%     subject to
%         dT == G_matrix * dS_tv(:)
% cvx_end
% toc;
%% save model
V1D = V_collec(:,:,end);
save([Tomo_path, '/Models/cross_well_att_1Dvp.mat'],'V1D')
%% plot
figure;
pcolor(Vtrue);
shading interp;
colormap(jet);
colorbar;
caxis([500 1100]);
axis ij;

figure;
pcolor(Vsm);
shading interp;
colormap(jet);
colorbar;
caxis([500 1100]);
axis ij;

figure;
pcolor(V_collec(:, :, end));
shading interp;
colormap(jet);
colorbar;
caxis([500 1100]);
axis ij;

figure;
pcolor(V_collec(:, :, end));
shading interp;
colormap(jet);
colorbar;
caxis([500 1100]);
axis ij;
hold on;
for i = 1:size(Raypoints, 1)
    plot(Raypoints{i}(:, 2), Raypoints{i}(:, 1), 'k');
    hold on;
end

figure;
pcolor(dS);
shading interp;
colormap(jet);
colorbar;
axis ij;

% figure;
% pcolor(dS_tv);
% shading interp;
% colormap(jet);
% colorbar;
% axis ij;

figure;
plot(Residual);

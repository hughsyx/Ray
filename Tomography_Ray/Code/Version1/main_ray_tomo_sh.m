%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Two Dimentional Tomography Using Matlab Functions
%                               Dongzhuo Li, SWP, Stanford
%                                           May 22, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
tic;

%% set paths
Tomo_path = '/Users/Dongzhuo/Documents/MATLAB/Tomography_Ray';
RT_path = '/Users/Dongzhuo/Documents/MATLAB/Eikonal';

%% marmousi test
load([Tomo_path, '/Models/marmousi_vp.mat']);
load([Tomo_path, '/Models/marmousi_vp_smoothed.mat']);
Vp = marmousi_vp;
Vp_smoothed = marmousi_vp_smoothed;
clear('marmousi_vp', 'marmousi_vp_smoothed');
Vp = Vp(21:120, 201:300);
Vp = modelSmooth(Vp, 3);
Vp_smoothed = Vp_smoothed(21:120, 201:300);

%% set parameter files
nDim = 2;
[nz, nx] = size(Vp);
ny = 0;
deltax = 1;
deltay = 1;
deltaz = 1;
xmin = 1;
ymin = 1;
zmin = 1;
vecZSource = 5:2:95;
vecXSource = 5*ones(1,length(vecZSource));
vecZReceiver = 5:2:95;
vecXReceiver = 95*ones(1,length(vecZReceiver));

%% compute traveltime map
AllTMaps = zeros(nz, nx, length(vecZSource));
for inxS = 1:length(vecZSource)
    AllTMaps(:, :, inxS) = msfm(Vp, [vecZSource(inxS); vecXSource(inxS)], true, false);
end

%% raytracing
count = 1;
for inxS = 1:length(vecZSource)
    for inxR = 1:length(vecZReceiver)
        [G_matrix(count,:), Raypoints{count}] = raytracing(AllTMaps(:, :, inxS),1:100,0,...
            1:100,1:100,0, 1:100, vecXSource(inxS), 0, vecZSource(inxS), ...
            vecXReceiver(inxR), 0, vecZReceiver(inxR));
        count = count + 1;
    end
end

% count = 1;
% for inxS = 1:length(vecZSource)
%     for inxR = 1:length(vecZReceiver)
%         Raypoints{count} = shortestpath(AllTMaps(:, :, inxS), [vecZReceiver(inxR); ...
%             vecXReceiver(inxR)], [vecZSource(inxS); vecXSource(inxS)],deltax,'euler');
%         count = count + 1;
%     end
% end

%% plot the rays
figure;
pcolor(Vp); shading interp; colormap(jet); axis ij;
hold on;
for i = 1:size(Raypoints, 2)
    plot(Raypoints{i}(1, :), Raypoints{i}(3, :), 'k');
    hold on;
end
toc;

%% need to find dT...
% Vp_temp = Vp_smoothed;
% for iTer = 1:3
%     AllTMaps = zeros(nz, nx, length(vecZSource));
%     for inxS = 1:length(vecZSource)
%         AllTMaps(:, :, inxS) = msfm(Vp_smoothed, [vecZSource(inxS); vecXSource(inxS)], true, false);
%     end
%     
%     count = 1;
%     for inxS = 1:length(vecZSource)
%         for inxR = 1:length(vecZReceiver)
%             [G_matrix(count,:), Raypoints{count}] = raytracing(AllTMaps(:, :, inxS),1:100,0,...
%                 1:100,1:100,0, 1:100, vecXSource(inxS), 0, vecZSource(inxS), ...
%                 vecXReceiver(inxR), 0, vecZReceiver(inxR));
%             count = count + 1;
%         end
%     end
% end
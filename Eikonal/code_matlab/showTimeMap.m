clear;
close all;
%% 2D
figure(1);
timeMap2D = load('timeMap2D.txt');
timeMap2D = reshape(timeMap2D, 425, 263);
imagesc(timeMap2D);
colormap('jet');
%% 3D
figure(2);
timeMap3D = load('timeMap3D.txt');
timeMap3D = reshape(timeMap3D, 10, 10, 10);
imagesc(timeMap3D(:,:,4));
colormap('jet');

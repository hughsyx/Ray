%%% This m-file is to build the initial 1D velocity model from traveltime
function [Sstart,Vstart] = build0Dvel0(sou_location,rec_location,...
    ForwardX,ForwardZ,DAOtime_obs)
%clear all;close all;clc;
%foldername = '../Experiments/velModelOneCircle';
%addpath('../commonCodes');
%load ([foldername,'/vmodel.mat']); %true velocity model
%load ([foldername,'/eikonalDAOdata.mat']); %load the traveltime of vtrue computed by eikonal solver
%%%%%%%%%%%%%%%%%%%%%%%%% Build Initial Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% velocity assigned on each node, source/receiver on the node %%%%%%%
Nsrc=size(sou_location,1);
Nrec=size(rec_location,1);
nx=length(ForwardX); nz=length(ForwardZ);
S0=0;
for i=1:Nsrc
    for j=1:Nrec
        S0 = S0 + DAOtime_obs(j,i)/norm(sou_location(i,:)-rec_location(j,:));
    end
end
S0 = S0/Nsrc/Nrec;
Sstart = S0*ones(nz,nx);
Vstart=1./Sstart;

%% in traveltime matrix, every row is one source, after reshape, it's receiver gather.
%%%%%%%%%% Plot Innitial Model %%%%%%%%
%{
figure;
subplot(121)
imagesc(ForwardX,ForwardZ,Sstart);
H=colorbar;set(get(H,'Title'),'string','s/m ');
hold on;
plot(sou_location(:,1),sou_location(:,2),'r*','markersize',8);
hold on;
plot(rec_location(:,1),rec_location(:,2),'g<','markersize',8);
title('Initial Slowness','FontSize',18,'fontweight','b');
xlabel('Distance (m)','FontSize',18);
ylabel('Depth (m)','FontSize',18);
axis image;
subplot(122)
imagesc(ForwardX,ForwardZ,Vstart);
H=colorbar;set(get(H,'Title'),'string','m/s ');
hold on;
plot(sou_location(:,1),sou_location(:,2),'r*','markersize',8);
hold on;
plot(rec_location(:,1),rec_location(:,2),'g<','markersize',8);
title('Initial Velocity','FontSize',18,'fontweight','b');
xlabel('Distance (m)','FontSize',18);
ylabel('Depth (m)','FontSize',18);
axis image;
drawnow;
%}
end

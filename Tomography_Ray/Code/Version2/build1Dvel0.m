%%% This m-file is to build the initial 1D velocity model from traveltime
function [Sstart,Vstart] = build1Dvel0(sou_location,rec_location,isou_location,...
    ForwardX,ForwardZ,DAOtime_obs)
%clear all;close all;clc;
%foldername = '../Experiments/velModelOneCircle';
%addpath('../commonCodes');
%load ([foldername,'/vmodel.mat']); %true velocity model
%load ([foldername,'/eikonalDAOdata.mat']); %load the traveltime of vtrue computed by eikonal solver
%%%%%%%%%%%%%%%%%%%%%%%%% Build Initial Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% velocity assigned on each node, source/receiver on the node %%%%%%%
Nsrc=size(sou_location,1);zsrc=sou_location(:,end);dz=ForwardZ(2)-ForwardZ(1);

z1=min(zsrc):dz:max(zsrc);
nx=length(ForwardX);nz=length(ForwardZ);
Sstart0=zeros(Nsrc,nx);
tra=diag(DAOtime_obs);
for j=1:Nsrc
    Sstart0(j,:)=tra(j)/(rec_location(j,1)-sou_location(j,1))*ones(1,nx);
end
[xx0,zz0]=meshgrid(ForwardX,zsrc);
[xx,zz1]=meshgrid(ForwardX,z1);
Sstart1=interp2(xx0,zz0,Sstart0,xx,zz1);% do interpolation for z within source range
% expend the starting model beyond source range
Sstart=zeros(nz,nx);
for j=1:min(isou_location(:,2))-1
    Sstart(j,:)=Sstart1(1,:);
end
for j=min(isou_location(:,2)):max(isou_location(:,2))
    Sstart(j,:)=Sstart1(j-min(isou_location(:,2))+1,:);
end
for j=max(isou_location(:,2))+1:nz
    Sstart(j,:)=Sstart1(end,:);
end
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

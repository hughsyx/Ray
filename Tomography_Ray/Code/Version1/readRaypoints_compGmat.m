function [T, Raypoints, G_matrix] = readRaypoints_compGmat(RT_path, nz, oz, dz, nx, ox, dx)
% read the computed raypoints file to obtain the travel time and the
% coordinates of those points on the corresponding every single ray
% Also compute the G matrix
% Dongzhuo Li @ Stanford
% May 21, 2016

rayID = fopen([RT_path, '/Out/Raypoints.txt'], 'r');

DATA = textscan(rayID, '%f%f%f%f');

rayIndex = find(isnan(DATA{3}) ~= 1);

T = DATA{4}(rayIndex);

Raypoints = cell(length(rayIndex), 1);

for inxRay = 1:length(rayIndex)-1
    Raypoints{inxRay} = [DATA{2}(rayIndex(inxRay)+1 : rayIndex(inxRay+1)-1) ...
        DATA{1}(rayIndex(inxRay)+1 : rayIndex(inxRay+1)-1)];
end
Raypoints{inxRay+1} = [DATA{2}(rayIndex(inxRay+1)+1 : end) ...
        DATA{1}(rayIndex(inxRay+1)+1 : end)];

%G_matrix = zeros(length(rayIndex), nz * nx);
G_matrix = sparse(length(rayIndex), nz * nx);
% dh = sqrt((Raypoints{1}(2,1)-Raypoints{1}(1,1))^2 ...
%     + (Raypoints{1}(2,2)-Raypoints{1}(2,2))^2);
   
for inxRay = 1:length(rayIndex)
    G_matrix(inxRay,:) = G_matrix_compute(Raypoints{inxRay}, nz, oz, dz, nx, ox, dx);
end
    
end
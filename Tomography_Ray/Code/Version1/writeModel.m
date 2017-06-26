function writeModel(RT_path, model)
% write model to the model file in the eikonal & raytracing package
% Dongzhuo Li @ Stanford
% May 21, 2016

[nz, nx] = size(model);

modelfileID = fopen([RT_path, '/Par/model2D.txt'], 'w');

for ix = 1:nx
    for iz =1:nz
        fprintf(modelfileID, '%d %d %f\n', ix, iz, model(iz, ix));
    end
end
fclose(modelfileID);
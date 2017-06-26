function writeParfile(RT_path, nDim, nx, ny, nz, deltax, deltay, deltaz, xmin, ymin, zmin, ...
    zSou, xSou, zRec, xRec)
% write parameters to the parameter file in the eikonal & raytracing
% package
% Dongzhuo Li @ Stanford
% May 21, 2016

%% write the Par_file
parfileID = fopen([RT_path, '/Par/Par_file.txt'], 'w');

fprintf(parfileID,'# DIMENSION\n');
fprintf(parfileID, '%d\n', nDim);

fprintf(parfileID,'\n');
fprintf(parfileID, '# NX\n');
fprintf(parfileID, '%d\n', nx);

fprintf(parfileID,'\n');
fprintf(parfileID, '# NY\n');
fprintf(parfileID, '%d\n', ny);

fprintf(parfileID,'\n');
fprintf(parfileID, '# NZ\n');
fprintf(parfileID, '%d\n', nz);

fprintf(parfileID,'\n');
fprintf(parfileID, '# DELTAX\n');
fprintf(parfileID, '%f\n', deltax);

fprintf(parfileID,'\n');
fprintf(parfileID, '# DELTAY\n');
fprintf(parfileID, '%f\n', deltay);

fprintf(parfileID,'\n');
fprintf(parfileID, '# DELTAZ\n');
fprintf(parfileID, '%f\n', deltaz);

fprintf(parfileID,'\n');
fprintf(parfileID, '# XMIN\n');
fprintf(parfileID, '%f\n', xmin);

fprintf(parfileID,'\n');
fprintf(parfileID, '# YMIN\n');
fprintf(parfileID, '%f\n', ymin);

fprintf(parfileID,'\n');
fprintf(parfileID, '# ZMIN\n');
fprintf(parfileID, '%f\n', zmin);

fprintf(parfileID,'\n');
fprintf(parfileID, '# NX_INVERSION\n');
fprintf(parfileID, '%d\n', nx);

fprintf(parfileID,'\n');
fprintf(parfileID, '# NY_INVERSION\n');
fprintf(parfileID, '%d\n', ny);

fprintf(parfileID,'\n');
fprintf(parfileID, '# NZ_INVERSION\n');
fprintf(parfileID, '%d\n', nz);

fprintf(parfileID,'\n');
fprintf(parfileID, '# DELTAX_INVERSION\n');
fprintf(parfileID, '%f\n', deltax);


fprintf(parfileID,'\n');
fprintf(parfileID, '# DELTAY_INVERSION\n');
fprintf(parfileID, '%f\n', deltay);

fprintf(parfileID,'\n');
fprintf(parfileID, '# DELTAZ_INVERSION\n');
fprintf(parfileID, '%f\n', deltaz);

fprintf(parfileID,'\n');
fprintf(parfileID, '# XMIN_INVERSION\n');
fprintf(parfileID, '%f\n', xmin);

fprintf(parfileID,'\n');
fprintf(parfileID, '# YMIN_INVERSION\n');
fprintf(parfileID, '%f\n', xmin);

fprintf(parfileID,'\n');
fprintf(parfileID, '# YMIN_INVERSION\n');
fprintf(parfileID, '%f\n', xmin);

fprintf(parfileID,'\n');
fprintf(parfileID, '# Velocity file\n');
fprintf(parfileID, './Par/model2D.txt\n');

fprintf(parfileID,'\n');
fprintf(parfileID, '# Source parameter file\n');
fprintf(parfileID, './Par/source_par.txt\n');

fprintf(parfileID,'\n');
fprintf(parfileID, '# Receiver parameter file\n');
fprintf(parfileID, './Par/receiver_par.txt\n');

fclose(parfileID);

%% write the source parameter file
souParID = fopen([RT_path, '/Par/source_par.txt'], 'w');
for i = 1:length(zSou)
    fprintf(souParID, '%d    %f    %f\n', i, zSou(i), xSou(i));
end
fclose(souParID);

%% write the receiver parameter file
recParID = fopen([RT_path, '/Par/receiver_par.txt'], 'w');
for i = 1:length(zSou)
    for j = 1:length(zRec)
        fprintf(recParID, '%d    %d    %f    %f\n', i, j, zRec(j), xRec(j));
    end
end
fclose(recParID);


end
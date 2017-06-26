nz = 60;
ny = 50;
nx = 80;

A = 2000*ones(nz,nx);
A(41:61, 41:61) = 100;
F = fspecial('gaussian', [5 5], 5);
A = conv2(A, F, 'same');

% A = A(:);

fid1 = fopen('../Par/model2D.txt', 'wt');
for i = 1:nx
    for j =1:nz
        fprintf(fid1, '%f %f %f\n', j, i, A(j,i));
    end
end
fclose(fid1);




% fid1 = fopen('./Par/model2D.txt', 'wt');
% for i = 1:nz
%     for j=1:nx
%         fprintf(fid1, '%f %f %f\n', i, j, A(i, j));
%     end
% end
% fclose(fid1);

B = 2000*ones(nz, nx, ny);
fid2= fopen('../Par/model3D.txt', 'wt');
for i = 1:ny
    for j = 1:nx
        for k = 1:nz
            fprintf(fid2, '%f %f %f %f\n', k, j, i, B(k, j, i));
        end
    end
end
fclose(fid2);

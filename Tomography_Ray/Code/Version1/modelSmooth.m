function model_sm = modelSmooth(model, fSize)
G = fspecial('gaussian', [fSize, fSize], 3);
% G = 1/dim_operator^2*ones(dim_operator,dim_operator);
% G= 1/10/dim_operator*ones(10,dim_operator);
model4sm = padarray(model, [fSize, fSize], 'replicate', 'both');
model_conv = conv2(full(model4sm), G, 'same');
model_sm = model_conv(fSize+1:size(model_conv,1) - fSize, ...
    fSize+1:size(model_conv,2) - fSize);
end
function [im4] = opencomps(imf)
for i1 = 1:8
    im = imread(imf, i1);
    im2 = reshape(im,[],1);
    %im3 = int16(im2 .* 1000);
    %im4(:,i1) = single(im3) ./ 1000;
    im4(:,i1) = im2;
end
end
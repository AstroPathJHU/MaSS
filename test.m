function test(i1)
mycols = [0 0 1;
        1 .7529 .7961;
        0 1 0;
        1 1 0;
        1 0 0;
        0.91 0.41 0.17;
        0 1 1;
        1 0 1;
        1 1 1;
        0 0 0;];
wd = '\\bki04\Segmentation\IF_Membrane\inform_data\Component_Tiffs';
fn = dir([wd,'\*seg.tif']);
%
fn1 = fullfile(fn(i1).folder, fn(i1).name);
%
for i1 = 1:10
    im1 = imread(fn1,i1);
    im_size = size(im1);
    im1 = reshape(im1,[],1);
    im(:,i1) = im1 / max(im1);
end
%
im2 = 200 .* asinh(1.5 .* im) * mycols;
ii = zeros(im_size(1) * im_size(2),1);
%
for i1 = 1:2
    im1 = imread(fn1, 13+i1);
    im1 = reshape(im1,[],1);
    ii = ii + (im1 > 0);
end
%
ii = find(ii);
im2(ii,:) = repmat([0 (1 * 255) 0], length(ii),1);
im2 = uint8(im2);
im2 = reshape(im2,[im_size, 3]);
im2 = im2(600:900,1000:1300, :);
%
imshow(im2)
end
    
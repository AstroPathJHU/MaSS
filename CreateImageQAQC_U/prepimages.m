function [im_dapi, im_nodapi] = prepimages(im, c_map, im_size, scol, seg)
%
% create dapi images first
%
im_dapi = 180 * sinh(1.5 * im) * c_map;
im_dapi(seg,:) = repmat([scol 0 0], length(seg),1);
im_dapi = uint8(im_dapi);
im_dapi = reshape(im_dapi,[im_size, 3]);
%
% create no dapi images next
%
im = im(:,2:end);
c_map = c_map(2:end,:);
im_nodapi = 180 * sinh(1.5 * im) * c_map;
if ~isempty(seg)
    im_nodapi(seg,:) = repmat([scol 0 0], length(seg),1);
end
im_nodapi = uint8(im_nodapi);
im_nodapi = reshape(im_nodapi,[im_size, 3]);
%
end
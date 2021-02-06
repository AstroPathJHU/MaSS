function []  = getsegfiles(Markers, q, fold, nam)
p = q;
%
% find the rows which do not have a membrane
%
mrows = isnan(p.MeanMembraneDAPI);
%
p = p(~mrows,:);
%
% make global cellids for the image before segmentation correction
%
CellID = p.('CellID');
CellID = double(CellID);
CellNum = p.('CellNum');
CellNum = double(CellNum);
Phenotype = p.('Phenotype');
%
% get segmentation outlines for all alternative segmentations
%
im3 = [];
in3 = [];
%
trows = false(length(CellID),length(Markers.altseg));
%
for i1 = 1:length(Markers.altseg)
    markalt = Markers.altseg{i1};
    %
    % get folder and image names for altseg
    %
    fdname = [fold,'Phenotyped\',markalt];
    iname = fullfile(fdname,[nam,'_binary_seg_maps.tif']);
    %
    % get rows of altseg cells
    %
    trows(:,i1) = strcmp(Phenotype, markalt);
    %
    % get inForm cellids of altseg cells
    %
    tcellnums = CellNum(trows(:,i1));
    tcellids = CellID(trows(:,i1));
    %
    % convert images to a single column vector and change cell nums to
    % cellids
    %
    % Membrane
    %
    im = imread(iname,4);
    %
    im = reshape(im,[],1);
    [a,b] = ismember(im,tcellnums); 
    ii2 = b(a,:);
    ii2 = tcellids(ii2);
    imn = zeros(size(im));
    imn(a,:) = ii2;
    %
    im3(:, i1 + 1) = imn;
    %
    % Nucleus
    %
    in = imread(iname, 2);
    %
    in = reshape(in,[],1);
    [a,b] = ismember(in,tcellnums); 
    ii2 = b(a,:);
    ii2 = tcellids(ii2);
    inn = zeros(size(in));
    inn(a,:) = ii2;
    %
    in3(:,i1 + 1) = inn;
    %
end
%
%get filenames for 1ry seg images
%
fdname = [fold,'Phenotyped\',Markers.seg{1}];
%
iname = fullfile(fdname,[nam,'_binary_seg_maps.tif']);
%
% get cellids of 1ry seg cells
%
trowsall = sum(trows,2) > 0;
cellids = CellID(~trowsall,1);
cellnums = CellNum(~trowsall,1);
%
% Membrane
%
im = imread(iname,4);
%
im = reshape(im,[],1);
[a,b] = ismember(im,cellnums);
ii2 = b(a,:);
ii2 = cellids(ii2);
imn = zeros(size(im));
imn(a,:) = ii2;
%
im3(:, 1) = imn;
%
% Nucleus
%
in = imread(iname, 2);
%
in = reshape(in,[],1);
[a,b] = ismember(in,cellnums);
ii2 = b(a,:);
ii2 = cellids(ii2);
inn = zeros(size(in));
inn(a,:) = ii2;
%
in3(:,1) = inn;
%
% get tissue seg
%
tisseg = imread(iname,1);
%
% get component_data
%
cfd = [fold, 'Component_Tiffs'];
cnm = [nam,'_'];
%
iname = fullfile(cfd, [cnm,'component_data.tif']);
%
props = imfinfo(iname);
for i1 = 1:(length(props) - 1)
    im = imread(iname, i1);
    cim(:,:,i1) = im; 
end
%
% print the images
%
imsize = [props(1).Height, props(1).Width];
ds.ImageLength = props(1).Height;
ds.ImageWidth = props(1).Width;
ds.Photometric = Tiff.Photometric.MinIsBlack;
ds.BitsPerSample   = props(1).BitDepth;
ds.SamplesPerPixel = 1;
ds.SampleFormat = Tiff.SampleFormat.IEEEFP;
ds.RowsPerStrip    = length(props(1).StripByteCounts);
ds.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
ds.Software = 'MATLAB';
ds.Compression = Tiff.Compression.LZW;
%
iname = fullfile(cfd, [cnm,'component_data_w_seg.tif']);
%
ii = Tiff(iname,'w');
for i3 = 1:(length(props) - 1)
    d = cim(:,:,i3);
    ii.setTag(ds);
    write(ii,d);
    writeDirectory(ii)
end
%
% Tissue
%
ii.setTag(ds)
tisseg = single(tisseg);
write(ii,tisseg);
%
% Nucleus
%
for i3 = 1:(size(in3, 2))
    writeDirectory(ii)
    d = reshape(in3(:,i3), imsize);
    d = single(d);
    ii.setTag(ds);
    write(ii,d);
end
%
% Membrane
%
for i3 = 1:(size(im3,2))
    writeDirectory(ii)
    d = reshape(im3(:,i3), imsize);
    d = single(d);
    ii.setTag(ds);
    write(ii,d);
end
%
close(ii)
end

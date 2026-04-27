%% mkimageid function
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2018
%% --------------------------------------------------------------
%% Description
%%% creates variables for a single image
%% --------------------------------------------------------------
%%
function [q, imageida, mycol, imc, simage] =...
    mkimageid(charts, inum, wd, Markers, doseg)
%
% set image output properties
%
imageida.ds.Photometric = Tiff.Photometric.RGB;
imageida.ds.BitsPerSample   = 8;
imageida.ds.SamplesPerPixel = 3;
imageida.ds.SampleFormat = Tiff.SampleFormat.UInt;
imageida.ds.RowsPerStrip = 41;
imageida.ds.MaxSampleValue = 256;
imageida.ds.MinSampleValue = 0;
imageida.ds.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
imageida.ds.Software = 'MATLAB';
imageida.ds.ResolutionUnit = Tiff.ResolutionUnit.Inch;
imageida.ds.XResolution = 300;
imageida.ds.YResolution = 300;
imageida.ds.Compression = Tiff.Compression.LZW;
%
% get chart that correspond to inum
%
nc = [charts(inum).folder,'\',charts(inum).name];
q = load(nc);
q = q.fData;
q.fname = charts(inum);
q.fig.CellXPos = q.fig.CellXPos + 1;
q.fig.CellYPos = q.fig.CellYPos + 1;
%
% some image designations
%
imageida.wd = wd;
imageida.id = extractBefore(q.fname.name,'cleaned_phenotype_table.mat');
%
% write out Tables that comes from this image
%
writetable(q.fig,[wd,'\Phenotyped\Results\QA_QC\Tables_QA_QC\',...
    erase(q.fname.name,'.mat'),'.csv']);
%
% image input fname for segmentation images
%
sim{1} = [wd,'\Phenotyped\',Markers.seg{1},'\',imageida.id];
for i1 = 1: length(Markers.altseg)
    sim{i1+1} = [wd,'\Phenotyped\',Markers.altseg{i1},'\',imageida.id];
end
%
% image output fname for the full Marker images
%
imageida.outfull = [wd,...
    '\Phenotyped\Results\QA_QC\Phenotype\All_Markers\',imageida.id];
%
% image output fname for lineage markers
%
for i1 = 1:length(Markers.lin)
    imageida.outABlin{i1} = [wd,...
        '\Phenotyped\Results\QA_QC\Phenotype\',Markers.lin{i1},'\',imageida.id];
    imageida.outABcoex{i1} = [wd,'\Phenotyped\Results\QA_QC\Lin&Expr_Coex\',...
        Markers.lin{i1},'\',imageida.id];
end
%
% image output fname name for additional lineage markers (ie coexpression)
% image output fname for expression marker coexpression on lineage markers
%
for i2 = 1:length(Markers.add)
    imageida.outABlin{i1+1} = [wd,'\Phenotyped\Results\QA_QC\Phenotype\',...
        Markers.add{i2},'\',imageida.id];
    imageida.outABcoex{i1+1} = [wd,'\Phenotyped\Results\QA_QC\Lin&Expr_Coex\',...
        Markers.add{i2},'\',imageida.id];
    i1 = i1+1;
end
%
% image output fname for expression markers
%
for i1 = 1: length(Markers.expr)
    imageida.outABexpr{i1} = [wd,'\Phenotyped\Results\QA_QC\Phenotype\',...
        Markers.expr{i1},'\',imageida.id];    
end
ii = ismember(Markers.all, Markers.expr);
imageida.exprlayer = Markers.Opals(ii);
%
idx = find(Markers.nsegs > 1);
idx_count = length(imageida.outABexpr);
%
if idx
    for i1 = 1:length(idx)
        cidx = idx(i1);
        for i2 = 2:Markers.nsegs(cidx)
            idx_count = idx_count + 1;
            str = [wd,'\Phenotyped\Results\QA_QC\Phenotype\',...
                Markers.all{cidx},'_',num2str(i2)];
            if ~exist(str, 'dir')
                mkdir(str);
            end
            imageida.outABexpr{idx_count} = [str,'\',imageida.id];
            imageida.exprlayer = [imageida.exprlayer;Markers.Opals(cidx)];
        end
    end
end
%
% fname for the component_Tiff image
%
iname = [wd,'\Component_Tiffs\',...
    imageida.id,'component_data.tif'];
%
% read in all component images
%
props = imfinfo(iname);
imageida.size = [props(1).Height, props(1).Width];
%
imageida.ds.ImageLength = props(1).Height;
imageida.ds.ImageWidth = props(1).Width;
%
% inForm often exports a fixed spectral deck (DAPI + 480/520/.../780) even
% when MergeConfig omits unused opals (e.g. no 480 or no 540). In that case
% grayscale page count exceeds length(Markers.Opals)+1 and we map by Opal.
%
imc = load_component_imc(iname, props, Markers);
%
mycol.all = Markers.mycol.all;
%
% lineage marker colors only
%
lins = ismember(Markers.all,Markers.lin);
mycol.lin = mycol.all(2:end-1,:);
mycol.lin = mycol.lin(lins,:);
%
% expression marker colors only
%
expr = ismember(Markers.all,Markers.expr);
mycol.expr = mycol.all(2:end-1,:);
mycol.expr = mycol.expr(expr,:);
%
%%%segmentation images%%%
%
if doseg
    %
    % get rows from each alternative segmentation in the main table
    %
    trows = false(height(q.fig),length(Markers.altseg));
    for i1 = 1:length(Markers.altseg)
        trows(:,i1) = strcmp(q.fig.Phenotype,Markers.altseg{i1});
        cellnums = double(q.fig.CellNum(trows(:,i1)));
        %
        % read in alternative segmentations; this only works if there is 
        % tissue segmentation and nuclear segmentation in the 
        % binary_seg image; cytoplasm
        %
        s1 = imread([sim{i1 + 1},'binary_seg_maps.tif'], 4);
        %
        % set cell labels of segmentation image that are not 
        % in the main table to zero
        %
        s1(~ismember(double(s1),cellnums)) = 0;
        %
        s1 = reshape(s1,[],1);
        simage3{i1 + 1} = s1;
    end
    %
    % get every row for alternative segmentation in the main table
    %
    trowsall = sum(trows,2) > 0;
    %
    % read in primary segmentation image
    %
    s1 = imread([sim{1},'binary_seg_maps.tif'],4);
    %
    % get cellnums of primary segmentation data
    % (ie data not in any alt segs)
    %
    cellnums = double(q.fig.CellNum(~trowsall,:));
    %
    s1(~ismember(double(s1),cellnums))=0;
    s1 = reshape(s1,[],1);
    %
    simage3{1} = s1;
    %
    % read in tissue segmentation
    %
    % sum the images across the segmentations to create a single unique
    % segmentation
    %
    simage2 = [simage3{:}];
    %
    simage = sum(simage2,2);
    %
    simage(simage>0) = .5;
    %
    simage = reshape(double(simage), imageida.size);
else
    simage = zeros(imageida.size);
end
%
end
%% load_component_imc
% Load normalized component stack columns: column 1 = DAPI, then one column
% per Markers.Opals entry (MergeConfig order). Supports (1) strict
% one-to-one page order when TIFF grayscale count matches MergeConfig, or
% (2) PhenoCycler-style full deck DAPI + [480 520 540 570 620 650 690 780]
% when the export contains all spectral slots but MergeConfig lists a subset.
function imc = load_component_imc(iname, props, Markers)
%
STANDARD_OPALS_AFTER_DAPI = [480, 520, 540, 570, 620, 650, 690, 780];
n_standard_deck = 1 + numel(STANDARD_OPALS_AFTER_DAPI);
%
is_gray = strcmp({props.ColorType}, 'grayscale');
gray_idx = find(is_gray);
if isempty(gray_idx)
    error(['No grayscale pages in component_data.tif: ', iname]);
end
%
ref_h = props(gray_idx(1)).Height;
ref_w = props(gray_idx(1)).Width;
same_hw = arrayfun(@(k) props(k).Height == ref_h && props(k).Width == ref_w, gray_idx);
gray_idx = gray_idx(same_hw);
n_gray = numel(gray_idx);
%
expected_cols = 1 + length(Markers.Opals);
%
if n_gray == expected_cols
    imc = read_gray_pages_sequential(iname, props, gray_idx, expected_cols);
    return
end
%
if n_gray == n_standard_deck && all(ismember(Markers.Opals, STANDARD_OPALS_AFTER_DAPI))
    warning('MaSS:ComponentTiffDeck:Remap', ...
        ['component_data.tif has the full %d-plane spectral stack (DAPI + 8 Opals); ', ...
        'mapping MergeConfig opals %s onto standard deck [DAPI %s].'], ...
        n_gray, mat2str(Markers.Opals(:)'), mat2str(STANDARD_OPALS_AFTER_DAPI));
    imc = zeros(ref_h * ref_w, expected_cols);
    imc(:, 1) = normalize_plane(imread(iname, gray_idx(1)));
    for k = 1:length(Markers.Opals)
        op = Markers.Opals(k);
        slot = find(STANDARD_OPALS_AFTER_DAPI == op, 1);
        if isempty(slot)
            error(['Opal ', num2str(op), ...
                ' is not in the standard PhenoCycler deck; cannot map component_data.tif']);
        end
        tif_page = gray_idx(1 + slot);
        imc(:, k + 1) = normalize_plane(imread(iname, tif_page));
    end
    return
end
%
% Detailed failure: explain mismatch for debugging (e.g. missing 480/540
% in merge vs export, or extra non-spectral pages).
%
layer_lines = cell(n_gray, 1);
for u = 1:n_gray
    k = gray_idx(u);
    layer_lines{u} = sprintf( ...
        '  page %u: %ux%u %s', k, props(k).Height, props(k).Width, props(k).ColorType);
end
bad_opals = Markers.Opals(~ismember(Markers.Opals, STANDARD_OPALS_AFTER_DAPI));
error(['Component image layers do not match active markers.\n', ...
    '  Grayscale pages (full-field): ', num2str(n_gray), ...
    '\n  Expected matrix columns (DAPI + MergeConfig opals): ', num2str(expected_cols), ...
    '\n  MergeConfig Opals (numeric): ', mat2str(Markers.Opals(:)'), ...
    '\n  Opals not in standard deck [480 520 540 570 620 650 690 780]: ', mat2str(bad_opals(:)'), ...
    '\n  Full-deck remap applies only when grayscale count == ', num2str(n_standard_deck), ...
    ' and every MergeConfig opal is in the standard list.\n', ...
    '  Grayscale pages:\n', strjoin(layer_lines, '\n')]);
end
%
function vec = normalize_plane(im_plane)
vec = double(reshape(im_plane, [], 1));
m = max(vec);
if m > 0
    vec = vec ./ m;
end
end
%
function imc = read_gray_pages_sequential(iname, props, gray_idx, ncols)
ref_h = props(gray_idx(1)).Height;
ref_w = props(gray_idx(1)).Width;
imc = zeros(ref_h * ref_w, ncols);
for j = 1:ncols
    tif_page = gray_idx(j);
    imc(:, j) = normalize_plane(imread(iname, tif_page));
end
end

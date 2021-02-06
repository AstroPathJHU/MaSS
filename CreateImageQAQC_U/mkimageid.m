function [q, imageid, mycol, imc, simage] =...
    mkimageid(charts, inum, wd, Markers, doseg)
%
% set image output properties
%
imageid.ds.Photometric = Tiff.Photometric.RGB;
imageid.ds.BitsPerSample   = 8;
imageid.ds.SamplesPerPixel = 3;
imageid.ds.SampleFormat = Tiff.SampleFormat.UInt;
imageid.ds.RowsPerStrip = 41;
imageid.ds.MaxSampleValue = 256;
imageid.ds.MinSampleValue = 0;
imageid.ds.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
imageid.ds.Software = 'MATLAB';
imageid.ds.ResolutionUnit = Tiff.ResolutionUnit.Inch;
imageid.ds.XResolution = 300;
imageid.ds.YResolution = 300;
imageid.ds.Compression = Tiff.Compression.LZW;
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
imageid.wd = wd;
imageid.id = extractBefore(q.fname.name,'cleaned_phenotype_table.mat');
%
% write out Tables that comes from this image
%
writetable(q.fig,[wd,'\Phenotyped\Results\QA_QC\Tables_QA_QC\',...
    erase(q.fname.name,'.mat'),'.csv']);
%
% image input fname for segmentation images
%
sim{1} = [wd,'\Phenotyped\',Markers.seg{1},'\',imageid.id];
for i1 = 1: length(Markers.altseg)
    sim{i1+1} = [wd,'\Phenotyped\',Markers.altseg{i1},'\',imageid.id];
end
%
% image output fname for the full Marker images
%
imageid.outfull = [wd,...
    '\Phenotyped\Results\QA_QC\Phenotype\All_Markers\',imageid.id];
%
% image output fname for lineage markers
%
for i1 = 1:length(Markers.lin)
    imageid.outABlin{i1} = [wd,...
        '\Phenotyped\Results\QA_QC\Phenotype\',Markers.lin{i1},'\',imageid.id];
    imageid.outABcoex{i1} = [wd,'\Phenotyped\Results\QA_QC\Lin&Expr_Coex\',...
        Markers.lin{i1},'\',imageid.id];
end
%
% image output fname name for additional lineage markers (ie coexpression)
% image output fname for expression marker coexpression on lineage markers
%
for i2 = 1:length(Markers.add)
    imageid.outABlin{i1+1} = [wd,'\Phenotyped\Results\QA_QC\Phenotype\',...
        Markers.add{i2},'\',imageid.id];
    imageid.outABcoex{i1+1} = [wd,'\Phenotyped\Results\QA_QC\Lin&Expr_Coex\',...
        Markers.add{i2},'\',imageid.id];
    i1 = i1+1;
end
%
% image output fname for expression markers
%
for i1 = 1: length(Markers.expr)
    imageid.outABexpr{i1} = [wd,'\Phenotyped\Results\QA_QC\Phenotype\',...
        Markers.expr{i1},'\',imageid.id];    
end
ii = ismember(Markers.all, Markers.expr);
imageid.exprlayer = Markers.Opals(ii);
%
idx = find(Markers.nsegs > 1);
idx_count = length(imageid.outABexpr);
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
            imageid.outABexpr{idx_count} = [str,'\',imageid.id];
            imageid.exprlayer = [imageid.exprlayer;Markers.Opals(cidx)];
        end
    end
end
%
% fname for the component_Tiff image
%
iname = [wd,'\Component_Tiffs\',...
    imageid.id,'component_data.tif'];
%
% read in all component images
%
props = imfinfo(iname);
imageid.size = [props(1).Height, props(1).Width];
%
imageid.ds.ImageLength = props(1).Height;
imageid.ds.ImageWidth = props(1).Width;
%
for i2 = 1:8
    if strcmp(props(i2).ColorType, 'grayscale')
        im(:,1) = reshape(imread(iname,i2),[],1);
        imc(:,i2) =(im(:,1)./max(im(:,1)));
    end
end
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
    simage = reshape(double(simage), imageid.size);
else
    simage = zeros(imageid.size);
end
%
end

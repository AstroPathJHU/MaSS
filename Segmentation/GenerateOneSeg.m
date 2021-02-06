%% GenerateOneSeg for Vectra data given inform & merge tables output
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2018
%% --------------------------------------------------------------
%% Will create the segmentation maps for mergetbls.m output
%%% This function removes cells not in a segmentation map and relabels the CellNums
%%% from inform output and replaces them with the corresponding unique CellIDs
%% --------------------------------------------------------------
%% input: 
%%% Must input three variables to the code
%%%      1. wd -->> *\inform_data directory
%%%      2. fname -->> filename
%%%      3. CoordTable -->> Coordinate table with 1 row 3 columns: filename;ImageX;ImageY
%% Usage: 
%%% wd =  'W:\Clinical_Specimen\M41_1\inform_data'; 
%%% fname = 'M41_1_[33870,10031]_cleaned_phenotype_table.csv';
%%% CoordTable = 
%%% ---------------- filename	            ImageX	 ImageY
%%% ---------------- M41_1_[5499,17961].csv	 5499	 17961
%%%
%% output: For each image in the case it will generate the following output
%%% \inform_data\Component_Tiffs\M41_1_[33870,10031]_component_data_w_seg.tif
%%% This is a 13 layer image stack; the first 8 layers are from the orginal
%%% component_data; the last 5 correspond to -> 1. Tissue Seg; 2. NuclearSeg_1; 
%%% 3. NuclearSeg_2; 4. MembraneSeg_1; 5. MembraneSeg_2
%%--------------------------------------------------------------
%%
function[] = GenerateOneSeg(wd,fname, CoordTable)
% 
formatspec = strcat({' %d16'},{' %s '},repmat('%d16 ',[1,3]),{ ' %s '},...
    repmat('%f32 ',[1,52]),{' '},repmat('%f32 ',[1,4]));
formatspec = formatspec{1};
%
%
ds.ImageLength = 1004;
ds.ImageWidth = 1344;
ds.Photometric = Tiff.Photometric.MinIsBlack;
ds.BitsPerSample   = 32;
ds.SamplesPerPixel = 1;
ds.SampleFormat = Tiff.SampleFormat.IEEEFP;
ds.RowsPerStrip    = 63;
ds.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
ds.Software = 'MATLAB';
%
% read in specified table
%
fdtin = [wd,'\Phenotyped\Results\Tables','\',fname];
f = readtable(fdtin,'Format',formatspec,'Delimiter',',');
%
% reformat the table and file name
%
T = [repmat(CoordTable,height(f),1),f(:,{'CellID','CellNum',...
    'Phenotype','CellXPos','CellYPos'})];
fname = extractBefore(fname, '_cleaned_phenotype_table.csv');
%
% read in 1st segmenation map
%
imnm1 = [wd,'\Phenotyped\PDL1\',fname,'_binary_seg_maps.tif'];
%
s1 = imread(imnm1,3);
n1 = imread(imnm1,2);
n1 = reshape(n1,[],1);
s1 = reshape(s1,[],1);
%
% get nontumor cellids and cellnums from table T
% 
ii = table2array(T(~strcmp(T.Phenotype,'Tumor'),'CellNum'));
ii2 = table2array(T(~strcmp(T.Phenotype,'Tumor'),'CellID'));
%
ii = double(ii);
ii2 = double(ii2);
%
% remove all cellnums not part of 'nontumor' phenotypes
%
s1(~ismember(s1,ii)) = 0;
n1(~ismember(n1,ii)) = 0;
%
% renumber the cellnums in segmenation maps to cellid values
%
s3 = zeros(length(s1),1);
n3 = zeros(length(n1),1);
%
for i2 = 1:length(ii)
    v = ii(i2);
    s3(ismember(s1,v)) = ii2(i2);
    n3(ismember(n1,v)) = ii2(i2);
end
%
% read in 2st segmenation map
%
imnm2 = [wd,'\Phenotyped\Tumor\',fname,'_binary_seg_maps.tif'];
%
s2 = imread(imnm2,3);
n2 = imread(imnm2,2);
t2 = imread(imnm2,1);
s2 = reshape(s2,[],1);
n2 = reshape(n2,[],1);
%
% get tumor cellids and cellnums from table T
% 
ii = table2array(T(strcmp(T.Phenotype,'Tumor'),'CellNum'));
ii2 = table2array(T(strcmp(T.Phenotype,'Tumor'),'CellID'));
%
ii = double(ii);
ii2 = double(ii2);
%
% remove all cellnums not part of 'nontumor' phenotypes
%
s2(~ismember(s2,ii)) = 0;
n2(~ismember(n2,ii)) = 0;
%
% renumber the cellnums in segmenation maps to cellid values
%
s4 = zeros(length(s2),1);
n4 = zeros(length(n2),1);
%
for i2 = 1:length(ii)
    v = ii(i2);
    s4(ismember(s2,v)) = ii2(i2);
    n4(ismember(n2,v)) = ii2(i2);
end
%
% reformat the label matricies to the image sizes and data type for export
%
s5 = [single(s3),single(s4)];
n5 = [single(n3),single(n4)];
t2 = single(t2);
%
s5 = reshape(s5,1004,1344,2);
n5 = reshape(n5,1004,1344,2);
%
% export the images; starting with the segmentation alone
%{
pathtout = [wd,'\Component_Tiffs'];
iname = [pathtout,'\',fname,'_seg.tif'];
%
ii = Tiff(iname,'w');
ii.setTag(ds);
write(ii,t2(:,:,1));
writeDirectory(ii)
ii.setTag(ds);
write(ii,n5(:,:,1));
writeDirectory(ii)
ii.setTag(ds);
write(ii,n5(:,:,2));
writeDirectory(ii)
ii.setTag(ds);
write(ii,s5(:,:,1));
writeDirectory(ii)
ii.setTag(ds);
write(ii,s5(:,:,2));
close (ii)
%}
% export the table for segmentation
%{
iname = [pathtout,'\',fname,'.csv'];
T = T(:,{'filename','ImageX','ImageY','CellID','Phenotype',...
    'CellXPos','CellYPos'});
writetable(T, iname)
%}
pathtout = [wd,'\Component_Tiffs'];
iname = [pathtout,'\',fname,'_component_data_w_seg.tif'];
imin = [pathtout,'\',fname,'_component_data.tif'];
c = readim(imin);
%
% export the segmentation matrix as addition layers to the component_data
%
ii = Tiff(iname,'w');
for i3 = 1:8
    d = c(:,:,i3);
    ii.setTag(ds);
    write(ii,d);
    writeDirectory(ii)
end
ii.setTag(ds);
write(ii,t2(:,:,1));
writeDirectory(ii)
ii.setTag(ds)
write(ii,n5(:,:,1));
writeDirectory(ii)
ii.setTag(ds)
write(ii,n5(:,:,2));
writeDirectory(ii)
ii.setTag(ds)
write(ii,s5(:,:,1));
writeDirectory(ii)
ii.setTag(ds)
write(ii,s5(:,:,2));
close(ii)
end

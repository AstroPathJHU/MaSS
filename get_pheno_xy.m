function [xy] = get_pheno_xy(filnm,marker,wd1,layers)
%%-----------------------------------------------------------
%% load the csv file with the given marker and extract xy positions of each
%% cell
%%-----------------------------------------------------------
%
% read in table
%
warning('off','MATLAB:table:ModifiedAndSavedVarnames')
%
% create format specifier for each column in the table
%
formatspec = strcat(repmat('%s ',[1,4]),{' '},repmat('%f32 ',[1,11]),...
    { ' %s '},repmat('%f32 ',[1,5]),{' '},repmat('%f32 ',[1,5*layers]),...
    { ' %s '},repmat('%f32 ',[1,5]),{' '},repmat('%f32 ',[1,5*layers]),...
    { ' %s '},repmat('%f32 ',[1,5]),{' '},repmat('%f32 ',[1,5*layers]),...
     { ' %s '},repmat('%f32 ',[1,4]),{' '},repmat('%f32 ',[1,5*layers]),...
    {' '},repmat('%s ',[1,2]),{' '}, repmat('%f32 ',[1,4]),{' '}, ....
    repmat('%s ',[1,2]));
formatspec = formatspec{1};
%
T = readtable([wd1,'\',marker,'\',filnm,'.txt'],'Format',formatspec,...
    'Delimiter','\t','TreatAsEmpty',{' ','#N/A'});
vars = T.Properties.VariableNames;
ii = find(contains(vars,'micron'));
if ii
    units{1} = 'microns';
else
    units{1} = 'pixels';
end
%
filnm2 = extractBetween(filnm,'[',']');
filnm2 = strsplit(filnm2{1},',');
fx = str2double(filnm2{1});
fy = str2double(filnm2{2});
%
fold = extractBefore(wd1,'Phenotyped');
fold = [fold,'Component_Tiffs'];
iname = [fold,'\',replace(filnm,...
    'cell_seg_data','component_data.tif')];
imageinfo = imfinfo(iname);
W = imageinfo.Width;
H = imageinfo.Height;
scalea = 10^4 *(1/imageinfo(1).XResolution);
if strcmp(units{1},'pixels')
    scale = 1;
elseif strcmp(units{1},'microns')
    scale = 10^4 *(1/imageinfo(1).XResolution);
end
%
xy = T(:,{'CellID','CellXPosition','CellYPosition','Phenotype'});
%
if strcmp(units{1},'microns')
    fx = (fx - scalea*(W/2)); %microns
    fy = (fy - scalea*(H/2)); %microns
elseif strcmp(units{1},'pixels')
    fx = (1/scalea .* fx - (W/2)); %pixels
    fy = (1/scalea .* fy - (H/2)); %pixles
end
%
if find(xy.CellXPosition > W)
    %
    xy.CellXPosition = 1/scale .* (xy.CellXPosition - fx);
    xy.CellYPosition = 1/scale .* (xy.CellYPosition - fy);
    %
    xy.CellXPosition = round(xy.CellXPosition);
    xy.CellYPosition = round(xy.CellYPosition);
    %
    ii = xy.CellXPosition < 1;
    xy.CellXPosition(ii) = 1;
    %
    ii = xy.CellYPosition < 1;
    xy.CellYPosition(ii) = 1;
end
%
end

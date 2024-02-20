%% function: getphenofield; 
%% --------------------------------------------------------------
%% Created by: Alex Szalay
% Last Edit: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% separate the different phenotype for each ABxx stuctures so that they
% only contain the positive phenotypes; get .tmp to be used for the Other
% classification based on the lowest Opal in the first segmentation
% algorithm; add the file name as an extra field in the output structure
%% --------------------------------------------------------
%%
function [f, e_code] = getphenofield(C, Markers, units)
e_code = 0;
%
% get the X and Y resolution + image size for unit conversions 
%
fold = extractBefore(C.fname.folder,'Phenotyped');
fold = [fold,'\Component_Tiffs'];
iname = [fold,'\',extractBefore(C.fname.name,...
    "]_cell_seg"),']_component_data.tif'];
if isempty(iname)
    iname = [fold,'\',extractBefore(C.fname.name,...
        "]_CELL_SEG"),']_component_data.tif'];
end
%
try
    imageinfo = imfinfo(iname);
catch
    strings = split(iname, '\');
    e_code = 21;
    f = [];
    return
end
W = imageinfo.Width;
H = imageinfo.Height;
scalea = 10^4 *(1/imageinfo(1).XResolution); % UM/PX
%
% get only the postive cells from each phenotype
%
for i3 = 1:length(Markers.all)
    %
    mark = Markers.all{i3};
    mark1 = lower(mark);
    %
    Markers.lall{i3} = mark1;
    %
    dat = C.(mark);
    dat2 = dat(strcmpi(dat.Phenotype,mark),:);
    %
    if strcmp(mark,'Tumor') && isempty(dat2)
        markb = Markers.all_original{i3};
        dat2 = dat(strcmpi(dat.Phenotype,markb),:);
    end
    %
    if ~isempty(dat2)
        %
        dat2.Phenotype = repmat({mark},height(dat2),1);
        %
        if strcmp(units{i3},'microns')
            fx = (dat2.fx - scalea*(W/2)); %microns
            fy = (dat2.fy - scalea*(H/2)); %microns
            dat2.CellXPos = floor(1/scalea .* (dat2.CellXPos - fx));
            dat2.CellYPos = floor(1/scalea .* (dat2.CellYPos - fy));
        end
        %
    end
    f.(mark1) = dat2;
end
%
% create a set of others
%
dat2 = C.(Markers.seg{1});
[~,loc] = ismember(Markers.seg, Markers.all);
%
if ~isempty(dat2)
    %
    dat2.Phenotype = repmat({'Other'},height(dat2),1);
    %
    if strcmp(units{loc}, 'microns')
        fx = (dat2.fx - scalea*(W/2)); %microns
        fy = (dat2.fy - scalea*(H/2)); %microns
        dat2.CellXPos = floor(1/scalea .* (dat2.CellXPos - fx));
        dat2.CellYPos = floor(1/scalea .* (dat2.CellYPos - fy));
    end
    %
end
f.tmp = dat2;
%
f.fname = C.fname;
end

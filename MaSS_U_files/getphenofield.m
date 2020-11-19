function f = getphenofield(C, Markers, units)
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
imageinfo = imfinfo(iname);
W = imageinfo.Width;
H = imageinfo.Height;
scalea = 10^4 *(1/imageinfo(1).XResolution); % UM/PX
[~, loc] = ismember(Markers.lin, Markers.all);
unit_name = repmat({'pixels'}, length(Markers.all),1);
unit_name(loc) = units;
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
        if strcmp(unit_name{i3},'microns')
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
[~,loc] = ismember(Markers.seg, Markers.lin);
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
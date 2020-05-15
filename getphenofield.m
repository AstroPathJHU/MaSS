function f = getphenofield(C, Markers, units)
%
% get the X and Y resolution + image size for unit conversions 
%
fold = extractBefore(C.fname.folder,'Phenotyped');
fold = [fold,'\Component_Tiffs'];
iname = [fold,'\',extractBefore(C.fname.name,...
        'cell_seg_data.txt'),'component_data.tif'];
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
        if strcmp(units{1},'microns')
            fx = (dat2.fx - scalea*(W/2)); %microns
            fy = (dat2.fy - scalea*(H/2)); %microns
        elseif strcmp(units{1},'pixels')
            fx = (1/scalea .* dat2.fx - (W/2)); %pixels
            fy = (1/scalea .* dat2.fy - (H/2)); %pixles
        end
        %
        if find(dat2.CellXPos > W) 
            dat2.CellXPos = 1/scale .* (dat2.CellXPos - fx);
            ii = dat2.CellXPos < 1;
            %
            dat2.CellXPos(ii) = 1;
            %
            dat2.CellYPos = 1/scale .* (dat2.CellYPos - fy);
            ii = dat2.CellYPos < 1;
            %
            dat2.CellYPos(ii) = 1;
        end
        %
    end
    f.(mark1) = dat2;
end
%
% create a set of others
%
dat2 = C.(Markers.seg{1});
%
if ~isempty(dat2)
    %
    dat2.Phenotype = repmat({'Other'},height(dat2),1);
    %
    if strcmp(units{1},'microns')
        fx = (dat2.fx - scalea*(W/2)); %microns
        fy = (dat2.fy - scalea*(H/2)); %microns
    elseif strcmp(units{1},'pixels')
        fx = (1/scalea .* dat2.fx - (W/2)); %pixels
        fy = (1/scalea .* dat2.fy - (H/2)); %pixles
    end
    %
    if find(dat2.CellXPos > W)
        dat2.CellXPos = 1/scale .* (dat2.CellXPos - fx);
        ii = dat2.CellXPos < 1;
        %
        dat2.CellXPos(ii) = 1;
        %
        dat2.CellYPos = 1/scale .* (dat2.CellYPos - fy);
        ii = dat2.CellYPos < 1;
        %
        dat2.CellYPos(ii) = 1;
    end
    %
end
f.tmp = dat2;
%
f.fname = C.fname;
end

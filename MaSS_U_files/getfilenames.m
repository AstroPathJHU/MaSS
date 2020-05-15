function [filenames, e]  = getfilenames(wd,Markers,uc)
%
%delete old Results if it exists & create a new results table
%
e = [];
%
filenames = cell(length(Markers.all),1);
nm = cell(1,length(Markers.all));
%
if exist([wd,'\Results\Tables'],'dir')
    try
        delete([wd,'\Results\Tables\*'])
    catch
        e = ['Error in ', uc,': could not delete old ',...
            '"*\Results\Tables folder";',...
            ' a file may be open or permission may be denied'];
        disp(e)
        return
    end
else
    mkdir (wd,'Results\Tables')
end
%
if exist([wd,'\Results\tmp_ForFiguresTables'],'dir')
    try
        delete([wd,'\Results\tmp_ForFiguresTables\*'])
    catch
        e = ['Error in ', uc,': could not delete old ',...
            '"*\Results\tmp_ForFiguresTables" folder;',...
            ' a file may be open or permission may be denied'];
        disp(e)
        return
    end
else
    mkdir (wd,'Results\tmp_ForFiguresTables')
end
%
% import file names for all inform file names
%
for i1 = 1:length(Markers.all)
    m = [wd,'\',Markers.all{i1},'\*cell_seg_data.txt'];
    filenames{i1,1} = dir(m);
    nm{:,i1} = {filenames{i1,1}(:).name};
end
%
%check that all fields are in all inform output
%
for i1 = 2: length(Markers.all)
    a(:,i1-1) = ismember(nm{:,1},nm{:,i1});
end
[x,~] = size(a);
ii = zeros(x,1);
%
ii(sum(a,2) == (length(Markers.all) - 1) , 1) = 1;
ii = logical(ii);
%
filenames = filenames{1,1}(ii);
%
% check segmentation 
%
if ~isempty(filenames)
    nm = extractBefore(filenames(1).name,'cell_seg_data.txt');
    %
    % get 1ry segmentation and see if it has proper layers
    %
    wd1 = [wd,'\',Markers.seg{1},'\'];
    iname = [wd1,nm,'binary_seg_maps.tif'];
    props = imfinfo(iname);
    if length(props) < 4
        e = ['Error in ', uc,...
            ': check binary segmentation maps',...
            ' in ',Markers.seg{1}];
        disp(e)
        return
    end
    %
    % check 2ry segmentations to see if they have proper layers
    %
    if ~isempty(Markers.altseg)
        for i1 = 1:length(Markers.altseg)
            mark = Markers.altseg{i1};
            wd1 = [wd,'\',mark,'\'];
            iname = [wd1,nm,'binary_seg_maps.tif'];
            props = imfinfo(iname);
            if length(props) < 4
                e = ['Error in ', uc,...
                    ': check binary segmentation maps',...
                    ' in ',mark];
                disp(e)
                return
            end
        end
    end
end
%
end
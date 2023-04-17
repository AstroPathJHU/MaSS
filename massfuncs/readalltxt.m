%% function: readalltxt; read all table for a given fname
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description:
% function takes in a filename, a data structure with marker information,
% and a directory location
% it loads the *cell_seg_data.txt for each marker output with the filename
% in the file tree described below in the input section
%% --------------------------------------------------------
%%
function [C, units, e_code] = readalltxt(filename,Markers,wd)
%%-------------------------------------------------------------------------
%% get all the phenotypes, coordinates, intensties into a single table
%% 
%%-------------------------------------------------------------------------
e_code = 0;
%
% Create a string vector which contains desired variables from the inform
% output and new variable names
%
if isfield(Markers, 'Membrane')
    layers = length(Markers.Opals) + 3;
else 
    layers = length(Markers.Opals) + 2;
end
%
% read the text files in each folder to a cell vector (read lineage and
% expr separate?)
%
[v,units] = cellfun(@(x) readtxt(...
    filename, x, wd, layers), Markers.lin, 'Uni', 0);
[v2,units2] = cellfun(@(x) readtxt(...
    filename, x, wd, layers), Markers.expr, 'Uni', 0);
%
[~, loc] = ismember(Markers.lin, Markers.all);
unit_name = repmat({'pixels'}, length(Markers.all),1);
unit_name(loc) = units;
%
[~, loc] = ismember(Markers.expr, Markers.all);
unit_name(loc) = units2;
units = unit_name;
%
% read in tables for multiple segmentation
%
idx = find(Markers.nsegs > 1);
idx_count = 0;
v3 = {};
%
if idx
   for i1 = 1:length(idx)
        cidx = idx(i1);
        for i2 = 2:Markers.nsegs(cidx)
            idx_count = idx_count + 1;
            x = [Markers.all{cidx},'_',num2str(i2)];
            try 
                [v3{idx_count},units3] = readtxt(filename,x, wd, layers); %#ok<AGROW>
            catch
                e_code = 19;
                return
            end
            
            if ~strcmp(units3, units(cidx))
                e_code = 16;
                return
            end
        end
   end
end
%
% merge cell vectors
%
v = [v,v2,v3];
%
% get the lineage tables out of v and put them into struct C
%
for i3 = 1:length(Markers.lin)
    C.(Markers.lin{i3}) = v{i3};
end
%
% get expr markers out of v and put them into struct C, merging tables for 
% antibodies with multiple segmentations
%
ii = ismember(Markers.all,Markers.expr);
nsegs = Markers.nsegs(ii);
i4 = 1;
i6 = 0;
%
for i3 = length(Markers.lin)+1:length(Markers.all)
    if nsegs(i4) > 1
        tbl = v{i3};
        for i5 = 2:nsegs(i4)
            i6 = i6 + 1;
            v3_count = length(Markers.all) + i6;
            tbl = [tbl;v{v3_count}]; %#ok<AGROW>
        end
        C.(Markers.expr{i4}) = tbl;
    else
        C.(Markers.expr{i4}) = v{i3};
    end
    i4 = i4+1;
end
%
% add filenames to data struct
%
C.fname = filename;
%
end

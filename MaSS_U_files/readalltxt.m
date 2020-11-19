function [C, units] = readalltxt(filename,Markers,wd)
%%-------------------------------------------------------------------------
%% get all the phenotypes, coordinates, intensties into a single table
%% 
%%-------------------------------------------------------------------------
%
% Create a string vector which contains desired variables from the inform
% output and new variable names
%
vars.Opals = Markers.Opals;
layers = length(Markers.Opals) + 2;
vars.names = cell(1,length(Markers.all)*8);
vars.type = [repmat({'Mean'},1,4), repmat({'Total'},1,4)];
vars.comp = repmat({'Nucleus', 'Membrane','EntireCell','Cytoplasm'},1,2);
%
% create DAPI column names -- there are additional vectors a2 - a4 are
% created so that the program is capable of handling different inform
% outputs
%
a(:,1)= cellfun(@(x,y)strcat(x,'DAPI','_DAPI_',...
    y,'_NormalizedCounts_TotalWeighting_'),vars.comp,vars.type,'Uni',0);
a2(:,1)= a(:,1);
%
% Sometimes DAPI is named, sometimes it isn't resulting in different
% column names which we use these vectors
%
a3(:,1)= cellfun(@(x,y)strcat(x,'DAPI',...
    y,'_NormalizedCounts_TotalWeighting_'),vars.comp,vars.type,'Uni',0);
a4(:,1)= a3(:,1);
%
b(:,1)= cellfun(@(x,y)strcat(y,x,'DAPI'),vars.comp,vars.type,'Uni',0);
%
% column names for all other markers
%
for i3 = 1:length(vars.Opals)
    %
    % tumor antibody should have been labeled with 'Tumor' in inForm
    %
    a(:,i3+1) = cellfun(@(x,y)strcat(...
        x,Markers.all{i3},'_Opal',num2str(vars.Opals(i3)),'_',y,...
        '_NormalizedCounts_TotalWeighting_'),vars.comp,vars.type,'Uni',0);
    %
    % However, mistakes happen so we check with the Tumor
    % antibody designation from inForm
    %
    a2(:,i3+1) = cellfun(@(x,y)strcat(...
        x,Markers.all_original{i3},'_Opal',num2str(vars.Opals(i3)),'_',y,...
        '_NormalizedCounts_TotalWeighting_'),vars.comp,vars.type,'Uni',0);
    %
    % Do the same for alternate DAPI naming
    %
    a3(:,i3+1) = a(:,i3+1);
    a4(:,i3+1) = a2(:,i3+1);
    %
    b(:,i3+1) = cellfun(@(x,y)strcat(y,x,num2str(vars.Opals(i3))),...
        vars.comp,vars.type,'Uni',0);
end
%
% put the strings into a single vector
%
zz = a';
zz = vertcat(zz(:))';
vars.namesin = cellfun(@(x)shortenName(x),zz,'Uni',0);
%
zz = a2';
zz = vertcat(zz(:))';
vars.namesin_original = cellfun(@(x)shortenName(x),zz,'Uni',0);
%
zz = a3';
zz = vertcat(zz(:))';
vars.namesin_dapierror = cellfun(@(x)shortenName(x),zz,'Uni',0);
%
zz = a4';
zz = vertcat(zz(:))';
vars.namesin_dapierror_original = cellfun(@(x)shortenName(x),zz,'Uni',0);
%
% New Table Names
%
zz = b';
zz = vertcat(zz(:))';
vars.names = zz;
%
% read the text files in each folder to a cell vector (read lineage and
% expr separate?)
%
[v,units] = cellfun(@(x) readtxt(...
    filename, x, vars, wd, layers), Markers.lin, 'Uni', 0);
[v2,~] = cellfun(@(x) readtxt(...
    filename, x, vars, wd, layers), Markers.expr, 'Uni', 0);
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
            [v3{idx_count},~] = readtxt(filename,x, vars, wd, layers);
        end
   end
end
%
% merge cell vectors
%
v = [v,v2,v3];
%
% convert names to something more understandable
%
for i4 = 1:length(v)
    v{i4}.Properties.VariableNames = ...
        [{'SampleName','SlideID','fx','fy',...
        'CellID','Phenotype','CellXPos',...
        'CellYPos','EntireCellArea'},vars.names{:}];
end
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
            tbl = [tbl;v{v3_count}];
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
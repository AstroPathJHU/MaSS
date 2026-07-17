%% function: file loop; merge all multipass tables for a given specimen
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/09/2019
%% --------------------------------------------------------------
%% Description
% merge all multipass tables for a given specimen using proper error
% handling
%% --------------------------------------------------------
%%
function [errors, nfiles, Markers] = fileloop(...
    wd, sname, filenms, sum_filenames, Markers, logstring, imall)
%
tic
errors = cell(length(filenms), 1);
nfiles = 0; %#ok<NASGU>
%
% attempt to start/recover pool; if unavailable continue serially
%
max_workers = min(12, max(1, ceil(length(filenms)/2)));
[use_parallel, ~, par_msg] = startparpool(max_workers, 'MaSS');
if ~isempty(par_msg)
    mywritetolog(wd, sname, logstring, par_msg, 2, 'Tables');
end
%
% loop through the mergetbls function for each sample with error catching
%
seg_markers = [Markers.seg, Markers.altseg];
dep_markers = setdiff(Markers.all, seg_markers);
for i1 = 1:length(Markers.nsegs)
    if Markers.nsegs(i1) > 1
        for i2 = 2:Markers.nsegs(i1)
            dep_markers = [dep_markers, strjoin([Markers.all(i1), '_', num2str(i2)], '')];
        end
    end
end
marker_order = [seg_markers, dep_markers];
header = ['Image', marker_order, 'Clusters'];
writecell(header, [wd, '\Phenotyped\Results\cells.csv'],'WriteMode','overwrite')
output = cell(length(filenms), 1);
% fix output stuff plz
for i1 = 1:length(filenms)
    errors = cell(length(filenms), 1);
    [output{i1}, errors] = mergeloop(...
    filenms, sum_filenames, i1, errors, Markers, wd, 0, sname,...
    logstring, seg_markers, dep_markers);
end
%
% Check if matlab files were created for QAQC step
% If not run without Tumor check
%
figtabledir = [wd, '\Phenotyped\Results\tmp_ForFiguresTables'];
if length(dir(figtabledir)) < 3
    writecell(header, [wd, '\Phenotyped\Results\cells.csv'],'WriteMode','overwrite')
    output = cell(length(filenms), 1);
    if use_parallel
        parfor i1 = 1:length(filenms)
            errors = cell(length(filenms), 1);
            [output{i1}, errors] = mergeloop(...
            filenms, sum_filenames, i1, errors, Markers, wd, 1, sname,...
            logstring, seg_markers, dep_markers);
        end
    else
        for i1 = 1:length(filenms)
            errors = cell(length(filenms), 1);
            [output{i1}, errors] = mergeloop(...
            filenms, sum_filenames, i1, errors, Markers, wd, 1, sname,...
            logstring, seg_markers, dep_markers);
        end
    end
end
for i1 = 1:length(output)
    writecell(output{i1}, [wd, '\Phenotyped\Results\cells.csv'],'WriteMode','append')
end
%
% close parallel loop
%
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    delete(poolobj);
end
%
% get final result measures
%
filenms = dir([wd,...
    '\Phenotyped\Results\Tables\*_cleaned_phenotype_table.csv']);
nfiles = length(filenms);
%
end

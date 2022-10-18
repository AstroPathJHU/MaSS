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
    wd, sname, filenms, Markers, logstring, imall)
%
tic
errors = cell(length(filenms), 1);
nfiles = 0; %#ok<NASGU>
%
% start the parpool if it is not open;
% attempt to open with local at max cores, if that does not work attempt 
% to open with BG1 profile, otherwise parfor should open with default
%
if isempty(gcp('nocreate'))
    try
        numcores = feature('numcores');
        %
        if numcores > length(filenms) 
            numcores = ceil(length(filenms)/2);
        end
        %
        if numcores > 12
            numcores = 12; %#ok<NASGU>
        end
        evalc('parpool("local",numcores)');
    catch
    end
end
%
% loop through the mergetbls function for each sample with error catching
%
parfor i1 = 1:length(filenms)
    errors = cell(length(filenms), 1);
    [errors] = mergeloop(...
    filenms, i1, errors, Markers, wd, 0, sname, logstring);
end
%
% Check if matlab files were created for QAQC step
% If not run without Tumor check
%
figtabledir = [wd, '\Phenotyped\Results\tmp_ForFiguresTables'];
if length(dir(figtabledir)) < 3
    parfor i1 = 1:length(filenms)
        errors = cell(length(filenms), 1);
        [errors] = mergeloop(...
        filenms, i1, errors, Markers, wd, 1, sname, logstring);
    end
end
%
% close parallel loop
%
poolobj = gcp('nocreate');
delete(poolobj);
%
% get final result measures
%
filenms = dir([wd,...
    '\Phenotyped\Results\Tables\*_cleaned_phenotype_table.csv']);
nfiles = length(filenms);
%
end
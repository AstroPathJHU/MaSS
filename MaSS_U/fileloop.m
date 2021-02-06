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
    %
    % current filename
    %
    fname = filenms(i1);
    nm = extractBefore(fname.name,"_cell_seg");
    if isempty(nm)
        nm = extractBefore(fname.name,"_CELL_SEG");
    end
    log_name = nm;
    %
    % try each image through the merge tables function
    %
    try
        [fData, e_code] = mergetbls(fname, Markers, wd, imall);
        errors{i1} = 0;
        %
        if isempty(fData) && e_code == 0
            e_code = 14;
        end
        %
        if e_code ~= 0
            %
            err_handl(wd, sname, logstring, log_name, e_code);
            %
            nm = [nm,'_cleaned_phenotype_table.csv'];
            ftd = fullfile([wd,'\Phenotyped\Results\Tables'],nm);
            %
            delete(ftd)
            errors{i1} = 1;
            %
        end
        %
    catch
       err_handl(wd, sname, logstring, log_name, 14);
       errors{i1} = 1;  
    end
    %
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
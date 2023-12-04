%% function: file loop; merge all multipass tables for a given specimen
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/09/2019
%% --------------------------------------------------------------
%% Description
% merge all multipass tables for a given specimen using proper error
% handling
%% --------------------------------------------------------
%%
function [errors] = mergeloop(...
    filenms, sum_filenames, i1, errors, Markers, wd, imall, sname, logstring, seg_markers, dep_markers)
    %
    % current filename
    %
    fname = filenms(i1);
    sum_fname = sum_filenames(i1);
    nm = extractBefore(fname.name,"_cell_seg");
    if isempty(nm)
        nm = extractBefore(fname.name,"_CELL_SEG");
    end
    log_name = nm;
    %
    % try each image through the merge tables function
    %
    try
        [fData, e_code, err_msg] = mergetbls(fname, sum_fname, Markers, wd, imall, i1, seg_markers, dep_markers);
        errors{i1} = 0;
        %
        if isempty(fData) && e_code == 0
            e_code = 14;
        end
        %
        if e_code ~= 0
            %
            disp(e_code);
            err_handl(wd, sname, logstring, log_name, e_code, 'Tables', err_msg);
            errors{i1} = 1;
            %
        end
        %
    catch EM
        disp(EM);
       err_handl(wd, sname, logstring, log_name, 14, 'Tables', '');
       errors{i1} = 1;  
    end
    %
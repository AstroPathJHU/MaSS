%% imageloop
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2018
%% --------------------------------------------------------------
%% Description
%%%  Loop through all images with proper error handling
%% --------------------------------------------------------------
%%
%% imageloop
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2018
%% --------------------------------------------------------------
%% Description
%%%  Loop through all images with proper error handling
%% --------------------------------------------------------------
%%
function err_val = imageloop(wd, uc, logstring, Markers, charts, doseg, use_parallel)
%
err_val = 0;
e = cell(length(charts),1);
if nargin < 7
    use_parallel = false;
end
%
if use_parallel
    parfor i2 = 1:length(charts)
        log_name = extractBefore(charts(i2).name, '_cleaned');
        disp(charts(i2).name);
        try
            err_str = ['CreateQAQC ', log_name, ' started'];
            mywritetolog(wd, uc, logstring, err_str, 2, 'QA_QC');
            [s{i2}, imageid{i2}, mycol{i2}, imc, simage]...
                =  mkimageid(charts, i2, wd, Markers, doseg); %#ok<PFOUS,PFBNS>
            s{i2}.fig.CellXPos = round(s{i2}.fig.CellXPos);
            s{i2}.fig.CellYPos = round(s{i2}.fig.CellYPos);
            [s{i2},ima, imas] = mkphenim(s{i2}, Markers, mycol{i2},...
                imageid{i2}, imc, simage, doseg);
            mkindvphenim(...
                s{i2}, mycol{i2}, imageid{i2},...
                imc, simage, Markers, ima, imas);
            mkexprim(...
                mycol{i2}, imageid{i2}, imc, Markers, ima, imas, wd);
            err_str = ['CreateQAQC ', log_name, ' finished'];
            mywritetolog(wd, uc, logstring, err_str, 2, 'QA_QC');
            e{i2} = 0;
        catch E
            if contains(E.message, "expression segmentation")
                e{i2} = 1;
                err_str = ['ERROR: ', E.message, ' does not match primary segmentation'];
                mywritetolog(wd, uc, logstring, err_str, 2, 'QA_QC');
            else
                e{i2} = 1;
                err_str = ['ERROR: CreateQAQC ', log_name, ' failed'];
                mywritetolog(wd, uc, logstring, err_str, 2, 'QA_QC');
            end
        end
    end
else
    for i2 = 1:length(charts)
        log_name = extractBefore(charts(i2).name, '_cleaned');
        disp(charts(i2).name);
        try
            err_str = ['CreateQAQC ', log_name, ' started'];
            mywritetolog(wd, uc, logstring, err_str, 2, 'QA_QC');
            [s{i2}, imageid{i2}, mycol{i2}, imc, simage]...
                =  mkimageid(charts, i2, wd, Markers, doseg); %#ok<PFOUS,PFBNS>
            s{i2}.fig.CellXPos = round(s{i2}.fig.CellXPos);
            s{i2}.fig.CellYPos = round(s{i2}.fig.CellYPos);
            [s{i2},ima, imas] = mkphenim(s{i2}, Markers, mycol{i2},...
                imageid{i2}, imc, simage, doseg);
            mkindvphenim(...
                s{i2}, mycol{i2}, imageid{i2},...
                imc, simage, Markers, ima, imas);
            mkexprim(...
                mycol{i2}, imageid{i2}, imc, Markers, ima, imas, wd);
            err_str = ['CreateQAQC ', log_name, ' finished'];
            mywritetolog(wd, uc, logstring, err_str, 2, 'QA_QC');
            e{i2} = 0;
        catch E
            if contains(E.message, "expression segmentation")
                e{i2} = 1;
                err_str = ['ERROR: ', E.message, ' does not match primary segmentation'];
                mywritetolog(wd, uc, logstring, err_str, 2, 'QA_QC');
            else
                fprintf('Error message: %s\n', E.message);
                fprintf('File: %s\n', E.stack(1).file);
                fprintf('Function: %s\n', E.stack(1).name);
                fprintf('Line: %d\n', E.stack(1).line);
                e{i2} = 1;
                err_str = ['ERROR: CreateQAQC ', log_name, ' failed'];
                mywritetolog(wd, uc, logstring, err_str, 2, 'QA_QC');
            end
        end
    end
end
%
% make pie chart figure with heatmap
% in a separate loop because these figures are made via matlab graphs,
% which is not supported by parfor
%
%disp(str1{2});
%{
try
    %
    for i2 = 1:length(charts)
        %
        mkfigs(s{i2},Markers,imageid{i2}, mycol{i2});
        %
    end
    %
catch
    e{1} = ['Error in: ', uc, '; figure output failed for Image QA/QC'];
    disp(e{1})
    return
end
%}
%
% remove tmp_ForFiguresTables
%
% poolobj = gcp('nocreate');
% delete(poolobj);
%
if any(cell2mat(e))
    err_val = 14;
    return
else
    try
        rmdir([wd,'\Results\tmp_ForFiguresTables'],'s');
    catch
    end
end
%
end

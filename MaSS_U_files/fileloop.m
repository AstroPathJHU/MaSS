function [errors, errors2,tim,Markers] = fileloop(wd,sname,MergeConfig)
%
tim = cell(4,1);
%
tim{1} = datestr(now,'dd-mmm-yyyy HH:MM:SS');
%
tic
errors2 = [];
errors = cell(1,1);
Markers = [];
%
% get Markers structure
%
try
    Markers = createmarks(MergeConfig);
catch
    errors2 = ['Error in ',sname, ': check Batch ID files.'];
    disp(errors2);
    return
end
%
if strcmp(Markers.Tumor,...
        ' MaSS must have only one "Tumor" designation')
     errors2 = ['Error in ',sname,':', Markers.Tumor, '.'];
     disp(errors2);
     return
end
%
% get all filenames for a directory
%
try
    [filenms,errors2] =  getfilenames(wd,Markers, sname);
catch
    errors2 = ['Error in ',sname, ': check inForm output files.'];
    disp(errors2);
    return
end
if ~isempty(errors2)
    disp(errors2);
    return
end
%
tim{2} = length(filenms);
%
errors = cell(length(filenms));
%
% start the parpool if it is not open;
% attempt to open with local at max cores, if that does not work attempt 
% to open with BG1 profile, otherwise parfor should open with default
%
if isempty(gcp('nocreate'))
    try
        numcores = feature('numcores');
        if numcores > 12
            numcores = 12;
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
    %
    % try each image through the merge tables function
    %
    try
        [fData,~] = mergetbls(fname,Markers,wd);
        if isempty(fData)
            str = ['Error in ',fname.name,...
                ': check inForm Cell Analysis output.']
            %
            nm = extractBefore(fname.name,'cell_seg');
            nm = [nm,'cleaned_phenotype_table.csv'];
            ftd = fullfile([wd,'\Results\Tables'],nm);
            %
            delete(ftd)
            disp(str);
            errors{i1} = fname.name;
        end
    catch
       disp(['Error in ',fname.name,...
           ': check inForm Cell Analysis output.']);
       errors{i1} = fname.name;  
    end
    %
end
%
% close parallel loop
%
poolobj = gcp('nocreate');
delete(poolobj);
%
tim{3} = toc;
%
filenms = dir([wd,'\Results\Tables\*_cleaned_phenotype_table.csv']);
tim{4} = length(filenms);
%
end
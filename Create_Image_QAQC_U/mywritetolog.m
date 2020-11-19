%% function: mywritetolog
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 07/07/2020
%% --------------------------------------------------------------
%% Description
% create the log and output resultant error messages
%% --------------------------------------------------------------
%% input:
% err_val = exit code value indicating different errors
% loc = location in main code block of log file message
% wd = working directory of current specimen up to inform_data\Phenotyped
% tim = contains different file and time information
%%
function err_val = mywritetolog(wd, uc, logstring, err_str, locs, version)
%
tim = datestr(now,'yyyy-mm-dd HH:MM:SS');
logf = [wd,'\Phenotyped\Results\QA_QC\QA_QC.log'];
logp = [wd,'\Phenotyped\Results\QA_QC'];
err_val = 0;
%
if locs == 1
    %
    if exist(logp,'dir')
        try
            rmdir(logp,'s')
            mkdir(logp)
        catch
            err_val = 9;
            warning(['ERROR IN QA QC path:', wd,' ', uc]);
            return
        end
    else
        mkdir(logp)
    end
    %
    % create first line of file
    %
    str = [logstring, uc, ';CreateQAQC started-v',version,';', tim, '\r\n'];
    %
    fileID = fopen(logf,'wt');
    fprintf(fileID,str);
    fclose(fileID);
    %
end
%
% for error or warning messages write the message out in the correct format
%
if locs == 2
    %
    str = [logstring, uc, ';',err_str,';', tim, '\r\n'];
    %
    if isfile(logf)
        fileID = fopen(logf,'a');
    else 
        fileID = fopen(logf,'wt');
    end
    %
    fprintf(fileID,str);
    fclose(fileID);
    %
end
%
end
function log = createlog(errors, errors2, wd, tim, Markers)
%
% get non-empty cells, ie the names of images that had errors
%
errors = errors(~cellfun('isempty',errors));
%
% create file & folder
%
logfd = [wd,'\Results\Tables'];
if ~exist(logfd, 'dir')
    mkdir(logfd)
end
logf = [wd,'\Results\Tables\MaSSLog.txt'];
%
% create first line of file
%
tim1 = tim{1};
str = ['Merging a Single Sample for inForm Cell Analysis tables started - ',...
    tim1,'\r\n'];
tim1 = num2str(tim{2});
if ~isempty(tim1)
    str = [str,'     ',tim1,' inForm Cell Analysis *cell_seg_data.txt tables',...
        ' detected. \r\n'];
end
%
% add to list of any errors that may have occured
%
if ~isempty(errors2)
    %
    % if the code never made it passed Createmarks the error message is as
    % follows
    %
    str = [str, errors2,'\r\n'];
else
    if ~isempty(errors)
        str = [str,'There was/were an error(s) on the following image(s): \r\n'];
        for i1 = 1: length(errors)
            %
            % write out errors to the file
            %
            fname = errors{i1};
            str = [str,'     Problem processing image "',fname,...
                '": Error in fileloop("',fname,...
                '"): could not find inForm Cell Analysis output. \r\n'];
        end
    else
        %
        % if there were no errors output this message in line2 instead
        %
        tim1 = num2str(tim{3});
        tim2 = num2str(tim{4});
        %
        str = [str,'     Successfully merged ',tim2,' image data tables in ',...
            tim1, ' secs. \r\n'];
        %
        tim1 = dir([wd,'\Results\tmp_ForFiguresTables\*.mat']);
        tim1 = num2str(length(tim1));
        %
        if ~isempty(Markers.Tumor{1})
            tim2 = ['     Criteria for temporary tables: ',...
                'fields with 600 total cells and 60 tumor cells.'];
        else
              tim2 = ['     Criteria for temporary tables: ',...
                'fields with 600 total cells.'];
        end
        str = [str,tim1,...
            ' temporary tables for QA/QC figures printed. \r\n',tim2,' \r\n'];
    end
end
dt = datestr(now,'dd-mmm-yyyy HH:MM:SS');
str = [str,...
    'Merging a single sample for inForm Cell Analysis tables complete - ',...
    dt,'\r\n'];
%
% now write out the string
%
fileID = fopen(logf,'wt');
fprintf(fileID,str);
fclose(fileID);
%
log = str;
%
end
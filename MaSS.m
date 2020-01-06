%% function: MaSS (Merge a Single Sample)
%% --------------------------------------------------------------
%% Created by: Benjamin Green and Alex Szalay - Johns Hopkins - 01/09/2019
%% --------------------------------------------------------------
%% Description
% run cleanup protocol for multipass phenotyping of a single sample. The
% output .csv files for each image is placed into *\inform output\Results
%% --------------------------------------------------------------
%% input: 
% This code is reliant on a Batch file as described in the BatchID
% Structure file, in a Batch folder one level above the 'specimen' folder.
% ---------!!!Please review the BatchID Structure !!!------------ 
% --------------!!! file before continuing !!!-------------------
%
% wd = working directory of current specimen up to
% *\inform_data\Phenotyped
% the *\inform_data directory should also contain a *\Component_Tiffs folder 
%
% The *\Phenotyped folder should contain: all 6 ABs that are analyzed by
% inform in seperate folders, designated by their Target names
%
%                  EX. \inform_data\Phenotyped\CD8
%                  EX. \inform_data\Phenotyped\Tumor
%                  EX. \inform_data\Phenotyped\PDL1
%                  EX. \inform_data\Phenotyped\PD1
%                  EX. \inform_data\Phenotyped\FoxP3
%                  EX. \inform_data\Phenotyped\CD163
% 
% Every ABxx folder should have cell_seg_data.txt for each image to be
% analyzed
%
% For segmentations; ABxx with the lowest Opal for each segmentation 
% type should have the *_binary_seg_maps.tif images exported
%
% The binary segmentation maps should contain 3 or 4 layers depending on
% whether or not TissueClassification is included; Note: when saving the
% inform algorithm be sure to turn on all mask layers to be sure all
% segmentations are exported
%
% the *_component_data.tif for each image should be in a directory one
% level up from the *\inform_data folder, in a Component_Tiffs folder
%
% sname = the sample name of the slide
%
%% Usage:
% wd = 'W:\Clincial_Specimen\M1_1\inform_data\Phenotyped'
% sname = 'M1_1'
% BatchIDxls = 'W:\Clincial_Specimen\BatchID\BatchID_01.xlsx'
%% output:
% *\inform_data\Results\Tables\*_cleaned_phenotyped_table.csv files 
% *\inform_data\Results\tmp_ForFiguresTables\* _cleaned_phenotyped_table.mat files 
%% --------------------------------------------------------
%%
function [] = MaSS(wd,sname,BatchIDxls)
%
% run the file loop for the sample
%
[errors,errors2,tim,Markers] = fileloop(wd,sname,BatchIDxls);
%
% create a error log file for the sample
%
createlog(errors,errors2,wd,tim,Markers);
%
end
%% function: file loop; merge all multipass tables for a given specimen 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/09/2019
%% --------------------------------------------------------------
%% Description
% merge all multipass tables for a given specimen using proper error
% handling
%% --------------------------------------------------------
%%
function [errors, errors2,tim,Markers] = fileloop(wd,sname,BatchIDxls)
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
    Markers = createmarks(BatchIDxls);
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
    [filenms,errors2] =  getfilenames(wd,Markers);
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
%% function: createmarks; Create Markers data structure
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% function takes in a folder location and creates the Markers data
% structure
%% --------------------------------------------------------
%%
function[Markers] = createmarks(BatchIDxls)
%
warning('off','MATLAB:table:ModifiedAndSavedVarnames')
%
try
    BIDtbl = readtable(BatchIDxls);
    B = BIDtbl(:,{'Opal','Target',...
    'TargetType','CoexpressionStatus','SegmentationStatus',...
    'SegmentationHierarchy', 'ImageQA_QC', 'NumberofSegmentations'});
catch
    return
end
%
% remove DAPI row
%
dr = strcmp(B.Opal, 'DAPI');
B(dr,:) = [];
%
% start setting up Markers struct
%
Markers.Opals = B.Opal;
Markers.all = B.Target;
%
ii = strcmp('Tumor',B.ImageQA_QC);
Markers.all_original = Markers.all;
%
% change Tumor marker designation to 'Tumor'
%
if sum(ii) == 1
    Markers.all(ii) = {'Tumor'};
    Markers.Tumor{1} = 'Tumor';
elseif sum(ii) > 1
    Markers.Tumor =  ' MaSS must have only one "Tumor" designation';
    return
else
     Markers.Tumor{1} = '';
end
%
% get lineage and expression markers
%
LT = strcmp(B.TargetType,'Lineage');
Markers.lin = Markers.all(LT);
%
ET = strcmp(B.TargetType,'Expression');
Markers.expr = Markers.all(ET);
%
% get the markers with multiple segmentations, this will only be a
% capability on expression markers
%
nsegs = B.NumberofSegmentations;
nsegs = cellfun(@(x) str2double(x), nsegs, 'Uni',0);
nsegs = cell2mat(nsegs);
if find(nsegs(~ET) ~= 1)
     Markers.nsegs =  [' MaSS can only handle expression markers',...
         ' with multiple segmentations'];
    return
end
Markers.nsegs = nsegs;
%
% Set up segmentation status to define number of segmentations and which is
% the primary segmentation
%
SS = B.SegmentationStatus;
SS = cell2mat(SS);
SS = str2num(SS);
SS = SS(nsegs == 1);
mn = Markers.all(nsegs == 1);
%
% get number of different segmentations, remove markers with multiple
% segmentations from the contention
%
[~,y,~] = unique(SS);
ii = y(1);
%
Markers.seg = mn(ii);
%
Markers.altseg = cell(length(y)-1,1);
for i1 = 2:length(y)
    ii = y(i1);
    Markers.altseg(i1-1) = mn(ii);
end
%
% get coexpression status for lineage markers
%
CS = B.CoexpressionStatus(LT);
ii = ~strcmp(CS,'NA');
CS = CS(ii);
%
% track the corresponding target
%
TCS = Markers.lin(ii);
%
% get segmentation heirarchy
%
SH = B.SegmentationHierarchy;
Markers.SegHie = SH(LT);
%
% CS that is not NA in lineage markers; find which coexpressions are
% acceptable
%
Markers.add = [];
sego = [];
for i1 = 1:length(CS)
    %
    % get current target and opal
    %
    T = TCS{i1};
    ii = strcmp(T,Markers.all);
    o = Markers.Opals(ii);
    o = o{1};
    %
    % check them against rest of targets in coexpression
    %
    CStest = CS(~strcmp(TCS,T));
    TCStest = TCS(~strcmp(TCS,T));
    %
    for i2 = 1:length(CStest)
        o1 = CStest{i2};
        T1 = TCStest{i2};
        %
        % if the current target matches one of the targets in the rest 
        %
        if contains(o1,o)
            %
            % if the Markers.add is not empty; are both markers already
            % contained together
            %
            if ~isempty(Markers.add) && sum(contains(Markers.add,T)...
                & contains(Markers.add,T1))
                continue
            else
                track = length(Markers.add) + 1;
                Markers.add{track} = [T1,T];
                ii = strcmp(T1, Markers.lin);
                seg1 = Markers.SegHie(ii);
                ii = strcmp(T, Markers.lin);
                seg2 = Markers.SegHie(ii);
                %
                seg = max([str2num(seg1{1}),str2num(seg2{1})]);
                sego{track} = num2str(seg);
            end
        end
    end
end
%
Markers.SegHie = [Markers.SegHie;sego'];
%
% get coexpression status for expression markers
%
CS = B.CoexpressionStatus(ET);
for i1 = 1:length(CS)
    T = CS{i1};
    T = reshape(T,3,[])';
    [s,~] = size(T);
    x = arrayfun(@(x)contains(Markers.Opals,T(x,:)),1:s,'Uni',0);
    x = horzcat(x{:});
    Markers.Coex{i1} = sum(x,2);
end
%
%
%
Markers.Opals = cell2mat(Markers.Opals);
Markers.Opals = str2num(Markers.Opals);
%
Markers.all = Markers.all';
Markers.all_original = Markers.all_original';
Markers.lin = Markers.lin';
Markers.expr = Markers.expr';
Markers.nsegs = Markers.nsegs';
Markers.seg = Markers.seg';
Markers.altseg = Markers.altseg';
Markers.SegHie = Markers.SegHie';
end
%% function: getfilenames; Get names of all files in a directory
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/18/2019
%% --------------------------------------------------------------
%% Description
% function takes in a folder location and a Marker struct defined above
% and gets the names of all inform tables
%% --------------------------------------------------------------
%%
function [filenames, e]  = getfilenames(wd,Markers,uc)
%
%delete old Results if it exists & create a new results table
%
e = [];
%
filenames = cell(length(Markers.all),1);
nm = cell(1,length(Markers.all));
%
if exist([wd,'\Results\Tables'],'dir')
    try
        delete([wd,'\Results\Tables\*'])
    catch
        e = ['Error in ', uc,': could not delete old ',...
            '"*\Results\Tables folder";',...
            ' a file may be open or permission may be denied'];
        disp(e)
        return
    end
else
    mkdir (wd,'Results\Tables')
end
%
if exist([wd,'\Results\tmp_ForFiguresTables'],'dir')
    try
        delete([wd,'\Results\tmp_ForFiguresTables\*'])
    catch
        e = ['Error in ', uc,': could not delete old ',...
            '"*\Results\tmp_ForFiguresTables" folder;',...
            ' a file may be open or permission may be denied'];
        disp(e)
        return
    end
else
    mkdir (wd,'Results\tmp_ForFiguresTables')
end
%
% import file names for all inform file names
%
for i1 = 1:length(Markers.all)
    m = [wd,'\',Markers.all{i1},'\*cell_seg_data.txt'];
    filenames{i1,1} = dir(m);
    nm{:,i1} = {filenames{i1,1}(:).name};
end
%
%check that all fields are in all inform output
%
for i1 = 2: length(Markers.all)
    a(:,i1-1) = ismember(nm{:,1},nm{:,i1});
end
[x,~] = size(a);
ii = zeros(x,1);
%
ii(sum(a,2) == (length(Markers.all) - 1) , 1) = 1;
ii = logical(ii);
%
filenames = filenames{1,1}(ii);
%
% check segmentation 
%
if ~isempty(filenames)
    nm = extractBefore(filenames(1).name,'cell_seg_data.txt');
    %
    % get 1ry segmentation and see if it has proper layers
    %
    wd1 = [wd,'\',Markers.seg{1},'\'];
    iname = [wd1,nm,'binary_seg_maps.tif'];
    props = imfinfo(iname);
    if length(props) < 4
        e = ['Error in ', uc,...
            ': check binary segmentation maps',...
            ' in ',Markers.seg{1}];
        disp(e)
        return
    end
    %
    % check 2ry segmentations to see if they have proper layers
    %
    if ~isempty(Markers.altseg)
        for i1 = 1:length(Markers.altseg)
            mark = Markers.altseg{i1};
            wd1 = [wd,'\',mark,'\'];
            iname = [wd1,nm,'binary_seg_maps.tif'];
            props = imfinfo(iname);
            if length(props) < 4
                e = ['Error in ', uc,...
                    ': check binary segmentation maps',...
                    ' in ',mark];
                disp(e)
                return
            end
        end
    end
end
%
end
%% function: createlog; Create error log for images processed in the fileloop
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/18/2019
%% --------------------------------------------------------------
%% Description
% create a log file in the *\Results directory that tracks any failed
% fields that were not mentioned in the inform_error logs
%% --------------------------------------------------------------
%%
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
%% function:mergetbls; merge tables function for inform output
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% function takes in a filename, a data structure with marker information,
% and a directory location
% it generates a *\Tables\ directory and merged inform files
%% --------------------------------------------------------------
%%
function [fData,logf] = mergetbls(fname,Markers,wd)
%
% read in data
%
logf ='read in data';
[C, p2, units] = readalltxt(fname,Markers,wd);
%
if p2 ~= length(Markers.all)
    fData = [];
    logf = 'path error in file';
    return
end
%
% select proper phenotypes
%
logf = [logf,';get phenotype fields'];
f = getphenofield(C, Markers, units);
%
% remove cells within X number of pixels in cells
%
logf = [logf,';find other cells'];
d = getdistinct(f,Markers);
%
% removes double calls in hierarchical style
%
logf = [logf,';resolve coexpressions of lineage markers'];
q = getcoex(d,Markers);
%
% get polygons from inform and remove other cells inside secondary 
% segmentation polygons
%
logf = [logf,';find segmentation, remove others in tumor'];
a = getseg(q,Markers);
%
% determine expression markers by cells that are in X radius
%
logf = [logf,';find expression marker expression'];
fData = getexprmark(a,Markers);
%
% save the data
%
logf = [logf,';save tables'];
ii = parsave(fData,fname,Markers,wd);
logf = [logf,';make table for this figure: ',ii];
%
end
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
function [C, p2, units] = readalltxt(filename,Markers,wd)
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
a3(:,1)= cellfun(@(x,y)strcat(x,'DAPI',...
    y,'_NormalizedCounts_TotalWeighting_'),vars.comp,vars.type,'Uni',0);
a4(:,1)= a3(:,1);
b(:,1)= cellfun(@(x,y)strcat(y,x,'DAPI'),vars.comp,vars.type,'Uni',0);
%
% column names for all other markers
%
for i3 = 1:length(vars.Opals)
    a(:,i3+1) = cellfun(@(x,y)strcat(...
        x,Markers.all{i3},'_Opal',num2str(vars.Opals(i3)),'_',y,...
        '_NormalizedCounts_TotalWeighting_'),vars.comp,vars.type,'Uni',0);
    %
    % Markers.all_original
    % check the orginial designation if first read fails with Tumor
    % replaced designation
    %
    a2(:,i3+1) = cellfun(@(x,y)strcat(...
        x,Markers.all_original{i3},'_Opal',num2str(vars.Opals(i3)),'_',y,...
        '_NormalizedCounts_TotalWeighting_'),vars.comp,vars.type,'Uni',0);
    %
    % to maintain compatability with inForm
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
vars.namesin_updatedinForm = cellfun(@(x)shortenName(x),zz,'Uni',0);
%
zz = a4';
zz = vertcat(zz(:))';
vars.namesin_updatedinForm_original = cellfun(@(x)shortenName(x),zz,'Uni',0);
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
[v,p,units] = cellfun(@(x) readtxt(filename,x,vars,wd,layers),Markers.lin,'Uni',0);
[v2,p2,~] = cellfun(@(x) readtxt(filename,x,vars,wd,layers),Markers.expr,'Uni',0);
%
% read in tables for multiple segmentation
%
idx = find(Markers.nsegs > 1);
idx_count = 0;
v3 = {};
p3 = {};
%
if idx
   for i1 = 1:length(idx)
        cidx = idx(i1);
        for i2 = 2:Markers.nsegs(cidx)
            idx_count = idx_count + 1;
            x = [Markers.all{cidx},'_',num2str(i2)];
            [v3{idx_count},p3{idx_count},~] = readtxt(filename,x, vars, wd, layers);
        end
   end
end
%
% check that all paths are the same
%
p = [p,p2];
p2 = extractAfter(p,'im3');
p2 = sum(ismember(p2,'\flatw'));
%
p2 = length(Markers.all);
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
%% function: readtxt; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description:
% reads in a single cell seg data file
%% --------------------------------------------------------
%%
function [Q,P, units] = readtxt(filnm,marker,vars,wd,layers)
%%-----------------------------------------------------------
%% load the csv file with the given marker
%%-----------------------------------------------------------
%
% read in table
%
warning('off','MATLAB:table:ModifiedAndSavedVarnames')
%
% create format specifier for each column in the table
%
formatspec = strcat(repmat('%s ',[1,4]),{' '},repmat('%f32 ',[1,11]),...
    { ' %s '},repmat('%f32 ',[1,5]),{' '},repmat('%f32 ',[1,5*layers]),...
    { ' %s '},repmat('%f32 ',[1,5]),{' '},repmat('%f32 ',[1,5*layers]),...
    { ' %s '},repmat('%f32 ',[1,5]),{' '},repmat('%f32 ',[1,5*layers]),...
     { ' %s '},repmat('%f32 ',[1,4]),{' '},repmat('%f32 ',[1,5*layers]),...
    {' '},repmat('%s ',[1,2]),{' '}, repmat('%f32 ',[1,4]),{' '}, ....
    repmat('%s ',[1,2]));
formatspec = formatspec{1};
%
T = readtable([wd,'\',marker,'\',filnm.name],'Format',formatspec,...
    'Delimiter','\t','TreatAsEmpty',{' ','#N/A'});
%
% Find Coordinates
%
f = strsplit(filnm.name,{'_[',']_cell_seg_data.txt'});
f = strsplit(f{2},',');
T.fx = repmat(str2double(f{1}),height(T),1);
T.fy = repmat(str2double(f{2}),height(T),1);
%
% Selects variables of interest from originial data table by name
%
P = T(1,'Path');
P = table2array(P);
P = P{1};
unitstype = [repmat({'pixels'},4,1),repmat({'microns'},4,1)];
%
for i1 = [1,5]
    %
    if i1 == 1
        areavariable = 'EntireCellArea_pixels_';
    else
        areavariable = 'EntireCellArea_squareMicrons_';
    end
    %
    namestry(i1,:) = [{'SampleName','SlideID','fx','fy',...
        'CellID','Phenotype','CellXPosition',...
        'CellYPosition'},areavariable,...
        vars.namesin{:}];
    namestry(i1+1,:) = [{'SampleName','SlideID','fx','fy',...
        'CellID','Phenotype','CellXPosition',...
        'CellYPosition'},areavariable,...
        vars.namesin_original{:}];
    namestry(i1+2,:) =  [{'SampleName','SlideID','fx','fy',...
        'CellID','Phenotype','CellXPosition',...
        'CellYPosition'},areavariable,...
        vars.namesin_updatedinForm{:}];
    namestry(i1+3,:) =  [{'SampleName','SlideID','fx','fy',...
        'CellID','Phenotype','CellXPosition',...
        'CellYPosition'},areavariable,...
        vars.namesin_updatedinForm_original{:}];
end
%
fail = 1;
icount = 1;
while fail == 1
    try
        Q = T(:,namestry(icount,:));
        units = unitstype{icount};
        fail = 0;
    catch
        fail = 1;
    end
    %
    icount = icount + 1;
    if icount > length(namestry(:,1))
        fail = 0;
    end
end
%
end
%% function:shortenName; shorten the names of the inForm columns
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% If the names of the columns in inForm are longer than 63 char then
% matlab will only read the first 63 char of the header and the search
% strings need to be adjusted as such. This function takes in a string,
% checks if it is longer than 63 chars and, if it is, cuts the string at
% 63 chars
%% --------------------------------------------------------------
%%
function out = shortenName(n)
out = n;
if length(out) > 63
    out = extractBefore(out,64);
end
end
%% function: getphenofield; 
%% --------------------------------------------------------------
%% Created by: Alex Szalay
% Last Edit: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% separate the different phenotype for each ABxx stuctures so that they
% only contain the positive phenotypes; get .tmp to be used for the Other
% classification based on the lowest Opal in the first segmentation
% algorithm; add the file name as an extra field in the output structure
%% --------------------------------------------------------
%%
function f = getphenofield(C, Markers, units)
%
% get the X and Y resolution + image size for unit conversions 
%
fold = extractBefore(C.fname.folder,'Phenotyped');
fold = [fold,'\Component_Tiffs'];
iname = [fold,'\',extractBefore(C.fname.name,...
    'cell_seg_data.txt'),'component_data.tif'];
imageinfo = imfinfo(iname);
W = imageinfo.Width;
H = imageinfo.Height;
scalea = 10^4 *(1/imageinfo(1).XResolution);
if strcmp(units{1},'pixels')
    scale = 1;
elseif strcmp(units{1},'microns')
    scale = 10^4 *(1/imageinfo(1).XResolution);
end
%
% get only the postive cells from each phenotype
%
for i3 = 1:length(Markers.all)
    %
    mark = Markers.all{i3};
    mark1 = lower(mark);
    %
    Markers.lall{i3} = mark1;
    %
    dat = C.(mark);
    dat2 = dat(strcmpi(dat.Phenotype,mark),:);
    %
    if strcmp(mark,'Tumor') && isempty(dat2)
        markb = Markers.all_original{i3};
        dat2 = dat(strcmpi(dat.Phenotype,markb),:);
    end
    %
    if ~isempty(dat2)
        %
        dat2.Phenotype = repmat({mark},height(dat2),1);
        %
        if strcmp(units{1},'microns')
            fx = (dat2.fx - scalea*(W/2)); %microns
            fy = (dat2.fy - scalea*(H/2)); %microns
        elseif strcmp(units{1},'pixels')
            fx = (1/scalea .* dat2.fx - (W/2)); %pixels
            fy = (1/scalea .* dat2.fy - (H/2)); %pixles
        end
        %
        if find(dat2.CellXPos > W) 
            dat2.CellXPos = 1/scale .* (dat2.CellXPos - fx);
            ii = dat2.CellXPos < 1;
            %
            dat2.CellXPos(ii) = 1;
            %
            dat2.CellYPos = 1/scale .* (dat2.CellYPos - fy);
            ii = dat2.CellYPos < 1;
            %
            dat2.CellYPos(ii) = 1;
        end
        %
    end
    f.(mark1) = dat2;
end
%
% create a set of others
%
dat2 = C.(Markers.seg{1});
%
if ~isempty(dat2)
    %
    dat2.Phenotype = repmat({'Other'},height(dat2),1);
    %
    if strcmp(units{1},'microns')
        fx = (dat2.fx - scalea*(W/2)); %microns
        fy = (dat2.fy - scalea*(H/2)); %microns
    elseif strcmp(units{1},'pixels')
        fx = (1/scalea .* dat2.fx - (W/2)); %pixels
        fy = (1/scalea .* dat2.fy - (H/2)); %pixles
    end
    %
    if find(dat2.CellXPos > W)
        dat2.CellXPos = 1/scale .* (dat2.CellXPos - fx);
        ii = dat2.CellXPos < 1;
        %
        dat2.CellXPos(ii) = 1;
        %
        dat2.CellYPos = 1/scale .* (dat2.CellYPos - fy);
        ii = dat2.CellYPos < 1;
        %
        dat2.CellYPos(ii) = 1;
    end
    %
end
f.tmp = dat2;
%
f.fname = C.fname;
end
%% function: getdistinct;
%% --------------------------------------------------------------
%% Created by: Alex Szalay
%%% Last Edit: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% define the Other cells;start by getting all cells from the lowest 
% expression marker opal in the 1ry segmentation
% eliminate all the cells within 6 pixels of a lineage classification,
% for 'Other' cells that are within 6 pixels of each other remove cells
% based on
% 1. whether the cell is to close to two cells 
% 2. based on second expression marker intensity
%% --------------------------------------------------------
%%
function q = getdistinct(d,Markers)
%
% check for cells within 6 pixels of 'Other' cells
%
w = 6;
%
d.tmp = [d.tmp; d.(lower(Markers.seg{1}))];
e = cellfun(@(x)min(edist(d.(lower(x)),d.tmp),[],1)',Markers.lin, 'Uni', 0);
%
for i1 = 1:length(Markers.lin)
    if isempty(e{i1})
        e{i1} = repmat(w+1,[height(d.tmp) 1]);
    end
end
%
e = [e{:}];
ii = sum(e <= w, 2);
ii = ii ~= 0;
d.tmp = d.tmp(~ii,:);
%
%
% remove the Others within 6 pixels of each other
%
r = edist(d.tmp,d.tmp);
r0 = (r<=6.0);
n0 = size(r0,1);
%
r0 = r0.*(1-eye(n0));
r0 = triu(r0);
[ix,iy] = find(r0>0);
%
% delete cells that are in the table more than once (ie a cell that is too
% close to two cells)
%
[~,b] = unique(iy);
rows = ismember(iy,iy(setdiff(1:length(iy),b)));
qx = unique(iy(rows));
d.tmp.Phenotype(qx,:) = {'Coll'};
iy = iy(~rows);
ix = ix(~rows);
%
% delete the rest of the cells based on the sum of expression marker profile
%
expr = Markers.Opals(ismember(Markers.all,Markers.expr));
ixe =  zeros(length(ix), length(Markers.expr));
iye =  zeros(length(iy), length(Markers.expr));
%
for i1 = 1:length(Markers.expr)
    opnm = num2str(expr(i1));
    tnm = ['MeanEntireCell',opnm];
    ixe(:,i1) = table2array(d.tmp(ix,tnm));
    iye(:,i1) =  table2array(d.tmp(iy,tnm));
end
%
rows =  sum(ixe,2) >= sum(iye,2);
%
qx = iy(rows);
d.tmp.Phenotype(qx,:) = {'Coll'};
%
qx = ix(~rows);
d.tmp.Phenotype(qx,:) = {'Coll'};
%
qx = ~strcmp(d.tmp.Phenotype,'Coll');
d.tmp = d.tmp(qx,:);
%
d.tmp.Phenotype = repmat(cellstr('Other'),height(d.tmp),1);
%
% make the next struct for coex%%%
%
q.other = d.tmp;
q.all = [];
for i3 = 1:length(Markers.lin)
    q.all = [q.all;d.(lower(Markers.lin{i3}))];
end
%
% keep original phenotype
%
q.og = [q.all; q.other];
%
q.fname = d.fname;
for i3 = 1:length(Markers.expr)
    mark = lower(Markers.expr{i3});
    q.(mark) = d.(mark);
end

end
%% function: getcoex; 
%% --------------------------------------------------------------
%% Created by: Alex Szalay
% Last Edit: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% resolve the cells with multiple lineage marker calls within 6 pixels of
% each other. This is done in three steps:
% 1. remove cells that are of the same lineage within 6 pixels of each 
% other based off of expression marker intensity;
% 2. second change accepted coexpression to [1st Opal Target name, 2nd 
% Opal Target name]; 
% 3. remove cells based off of the segmentation hierarchy in the table 
%% --------------------------------------------------------------
%%
function q = getcoex(d,Markers)
%
% compute distance matrix with overlapping cells
%
[c] = createoverlapmatrix(d.all);
%
q = d;
q.c = c;
%
% remove missegmented cells
%
[q,c] = remvover(q,c,Markers);
%
% find coexpression after removal of missegmented cells & replace those
% calls in segmentation table 
%
[q,~] = addcoex(q,c,Markers);
%
% recompute the distance matrix
%
[c2] = createoverlapmatrix(q.all);
%
% patch for multiple coexpressions
%
[q,~,Markers] = addmultiplecoex(q,c2,Markers);
%
% recompute the distance matrix
%
[c2] = createoverlapmatrix(q.all);
%
% remove double calls coex
%
[q,~] = remvcoex(q,c2,Markers);
%
q.fig = [q.all;q.other];
qx = ~strcmp(q.fig.Phenotype,'Coll');
q.fig = q.fig(qx,:);
%
q = rmfield(q,'other');
%
% check again there are any cells too close in the pixel boundary
%
[c2] = createoverlapmatrix(q.fig);
%
if ~isempty(c2)
    error('Cells within limit, 6 pixels, exist after cleanup')
end
%
end
%% function: createoverlapmatrix;
%% --------------------------------------------------------------
%% Created by: Alex Szalay - Johns Hopkins 
%% Editied by: Benjamin Green Johns Hopkins - 11/15/2019
%% --------------------------------------------------------------
%% Description
% function to compute the cells within 6 pixels of each other
%% --------------------------------------------------------------
%%
function [c] = createoverlapmatrix(d)
%
% compute distance matrix
%
e  = edist(d,d);
%
% filter on e<=6pixels
%
e0 = (e<=6.0);
n02 = size(e0,1);
%
% remove diagonal, and get upper triangle
%
e0 = e0.*(1-eye(n02));
e0 = triu(e0);
[ix2,iy2] = find(e0>0);
%
% assemble it to coex
%
x = d.CellXPos;
y = d.CellYPos;
p = d.Phenotype;
%
c = table(x(ix2),y(ix2),p(ix2),x(iy2),y(iy2),p(iy2),ix2,iy2);
c.Properties.VariableNames = {'x1','y1','p1','x2','y2','p2','i1','i2'};
c.dist = sqrt(double((c.x1-c.x2).^2+(c.y1-c.y2).^2));
%
end
%% function: remvover;
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% remove cells of the same type which are too close together based off
% of:
% 1. whether the cell is a part of an acceptable coexpression
% 2. based off of total expression marker intensity 
%% -------------------------------------------------------------
%%
function[q,c] = remvover(q,c,Markers)
%
% find cells of same type which are too close together
%
rows = strcmp(c.p1,c.p2);
q.flags.misseg = c(rows,:);
q.flags.misseg.c = find(rows>0);
%
for i2 = 1:length(Markers.add)
ii = Markers.add(i2);
%
SW = cellfun(@(x)startsWith(ii,x), Markers.all,'Uni',0);
SW = [SW{:}];
SWm = Markers.all{SW};
%
EW = cellfun(@(x)endsWith(ii,x), Markers.all,'Uni',0);
EW = [EW{:}];
EWm = Markers.all{EW};
%
% find coexpression before removal of same type cells
%
rows = (strcmp(c.p1,SWm) & strcmp(c.p2,EWm))|...
    (strcmp(c.p1,EWm) & strcmp(c.p2,SWm));
q.coex = c(rows,:);
q.coex.c = find(rows>0);
%
midx1 = q.flags.misseg.i1;
cidx1 = q.coex.i1;
midx2 = q.flags.misseg.i2;
cidx2 = q.coex.i2;
%
rows1 = ismember(midx1, cidx1)|ismember(midx1,cidx2);
rows2 = ismember(midx2, cidx1)|ismember(midx2,cidx2);
%
% if cells are too close together keep the cell if it is coepressed with
% another cell and delete the one that is not coexpressed
%
rows = rows1==1 & rows2==0;
%
q.all.Phenotype(q.flags.misseg.i2(rows))= {'Coll'};
c.p2(q.flags.misseg.c(rows)) = {'Coll'};
q.flags.misseg.p2(rows) = {'Coll'};
%
rows =  rows1 == 0 & rows2 == 1;
q.all.Phenotype(q.flags.misseg.i1(rows))= {'Coll'};
c.p1(q.flags.misseg.c(rows)) = {'Coll'};
q.flags.misseg.p1(rows) = {'Coll'};
%

end

[~,ii] = ismember(Markers.expr,Markers.all);
exprO = Markers.Opals(ii);
%
for i1 = 1:length(Markers.all)
    mark = zeros(length(Markers.all),1);
    mark(i1) = 1;
    ii = zeros(length(Markers.all), length(Markers.expr));
    for i2 = 1:length(Markers.Coex)
        ii(:,i2) = mark & Markers.Coex{i2};
    end
    ii = sum(ii,1);
    ii = logical(ii);
    %
    if sum(ii)
        %
        expr = arrayfun(@(x)['MeanMembrane',num2str(x)],exprO,'Uni',0);
        %
        mark = Markers.all{i1};
        rows = find(strcmp(q.flags.misseg.p1,mark));
        %
        t = table2array(q.all(:,q.all.Properties.VariableNames(expr)));
        %
        % get i1 expression total
        %
        midx1 = q.flags.misseg.i1(rows);
        eidx1 = t(midx1,ii);
        eidx1 = sum(eidx1,2);
        %
        % get i2 expression total
        %
        midx2 = q.flags.misseg.i2(rows);
        eidx2 = t(midx2,ii);
        eidx2 = sum(eidx2,2);
        %
        % get rows where i1 expression is more than i2 expression
        %
        rows4 = eidx1 > eidx2;
        rows1 = rows(rows4);
        %
        q.all.Phenotype(q.flags.misseg.i2(rows1))= {'Coll'};
        c.p2(q.flags.misseg.c(rows1)) = {'Coll'};
        q.flags.misseg.p2(rows1) = {'Coll'};
        %
        %
        % check that missegmented cell is not duplicated/ causing other
        % invalid coexpression
        %
        ii1 = q.flags.misseg.c(rows1);
        r = ismember(c.i2,c.i2(ii1));
        c.p2(r) = {'Coll'};
        r = ismember(c.i1,c.i2(ii1));
        c.p1(r) = {'Coll'};
        %
        rows1 = rows(~rows4);
        %
        q.all.Phenotype(q.flags.misseg.i1(rows1))= {'Coll'};
        c.p1(q.flags.misseg.c(rows1)) = {'Coll'};
        q.flags.misseg.p1(rows1) = {'Coll'};
        %
        % check that missegmented cell is not duplicated/ causing other
        % invalid coexpression
        %
        ii1 = q.flags.misseg.c(rows1);
        r = ismember(c.i2,c.i1(ii1));
        c.p2(r) = {'Coll'};
        r = ismember(c.i1,c.i1(ii1));
        c.p1(r) = {'Coll'};
    end
end
end
%% function: addcoex; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% find the acceptable coexpression after cells that are too close
% together are resolved; rename the cells to what the acceptable
% coexpression name is in the table; remove the cell which the acceptable
% coexpression marker designation in Tables ends with and keep the cell
% which the acceptable coexpression marker designation in Tables starts
% with
%% --------------------------------------------------------------
%%
function [q,c] = addcoex(q,c,Markers)
q.coex = table();
for i2 = 1:length(Markers.add)
    ii = Markers.add(i2);
    %
    SW = cellfun(@(x)startsWith(ii,x), Markers.all,'Uni',0);
    SW = [SW{:}];
    SWm = Markers.all{SW};
    %
    EW = cellfun(@(x)endsWith(ii,x), Markers.all,'Uni',0);
    EW = [EW{:}];
    EWm = Markers.all{EW};
    %
    % find coexpression and replace whatever the additional call starts 
    % with(ie. lower Opal); remove whatever the additional call ends with
    % (ie. higher Opal)
    %
    rows = strcmp(c.p1,SWm) & strcmp(c.p2,EWm);
    if sum(rows) > 0
        q.all.Phenotype(c.i1(rows))= {'Coll'};
        q.all.Phenotype(c.i2(rows))= ii;
        q.coex = [q.coex;c(rows,:)];
        c.p1(rows) = {'Coll'};
        c.p2(rows) = ii;
    end
    %
    rows = strcmp(c.p1,EWm) & strcmp(c.p2,SWm);
    if sum(rows) > 0
        q.all.Phenotype(c.i2(rows))= {'Coll'};
        q.all.Phenotype(c.i1(rows))= ii;
        q.coex = [q.coex;c(rows,:)];
        c.p2(rows) = {'Coll'};
        c.p1(rows) = ii;
    end
    %
end
end
%% function: remvcoex; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% remove the cells within 6 pixels of each which have conflicting lineage
% marker designations (lineage markers not in acceptable coexpression)
% based off of the segmentation hierarchy
%% --------------------------------------------------------------
%%
function [q,c2] = remvcoex(q,c2,Markers)
%
seg = cellfun(@str2double,Markers.SegHie);
seg1 = seg;
%
allmarkers = [Markers.lin, Markers.add];
while ~isempty(seg1)
    %
    s1 = max(seg1);
    seg1(seg1 == s1) = [];
    %
    ii = seg == s1;
    % 
    mark = allmarkers(ii);
    %
    for i1 = 1:length(mark)
        marka = mark{i1};
        rows = strcmp(c2.p2,marka) & ~strcmp(c2.p1,'Coll');
        q.all.Phenotype(c2.i2(rows))= {'Coll'};
        q.flags.(marka) = c2(rows,:);
        c2.p2(rows) = {'Coll'};
        %
        rows = strcmp(c2.p1,marka) & ~strcmp(c2.p2,'Coll');
        q.all.Phenotype(c2.i1(rows)) = {'Coll'};
        q.flags.(marka) = c2(rows,:);
        c2.p1(rows) = {'Coll'};
    end
    %
end
end
%% function: addmultiplecoex;
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 11/15/2019
%% --------------------------------------------------------------
%% Description
% for coexpression more than 2 cells 
%% --------------------------------------------------------------
%%
function [q,c2,Markers] = addmultiplecoex(q,c2,Markers)
%
seg = cellfun(@str2double,Markers.SegHie);
seg1 = seg;
%
allmarkers = [Markers.lin, Markers.add];
while ~isempty(seg1)
    %
    s1 = min(seg1);
    seg1(seg1 == s1) = [];
    %
    ii = seg == s1;
    % 
    mark = allmarkers(ii);
    if numel(mark) > 4
        ii = ~ismember(mark,Markers.add);
        markb = flip(mark(ii));
        markb = strcat(markb{:});
        mark = mark(~ii);
        %
        %
        for i1 = 1:length(mark)
            marka = mark{i1};
            rows = strcmp(c2.p2,marka) & ismember(c2.p1,mark);
            q.all.Phenotype(c2.i1(rows))= {'Coll'};
            q.all.Phenotype(c2.i2(rows))= {markb};
            q.flags.(marka) = c2(rows,:);
            c2.p1(rows) = {'Coll'};
            c2.p2(rows) = {markb};
            %
            rows = strcmp(c2.p1,marka) & ismember(c2.p2,mark);
            q.all.Phenotype(c2.i2(rows)) = {'Coll'};
            q.all.Phenotype(c2.i1(rows))= {markb};
            q.flags.(marka) = c2(rows,:);
            c2.p2(rows) = {'Coll'};
            c2.p1(rows) = {markb};
        end
        Markers.add = [Markers.add,markb];
        Markers.SegHie = [Markers.SegHie,num2str(s1)];
        %
    end
    %
end
end
%% function: getseg;
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% remove other cells with cell centers within non primary segmentation 
% boundaries; If a cell is on the boundary of the tissue the algorithm
% will close the cell to complete the calculation
%% --------------------------------------------------------------
%%
function [p] = getseg(q,Markers)
p = q;
%
% find the rows which do not have a membrane
%
mrows = isnan(p.fig.MeanMembraneDAPI);
nmr = p.fig(mrows,:);
p.fig = p.fig(~mrows,:);
%
% make global cellids for the image before segmentation correction
%
CellID = (1:1:height(p.fig))';
p.fig.Properties.VariableNames('CellID') = {'CellNum'};
p.fig = [table(CellID),p.fig];
%
% find the other cells before the tumor segmentation correction
%
orows = find(strcmp(p.fig.Phenotype, 'Other'));
%
% get segmentation outlines for all alternative segmentations
%
im3 = cell(length(Markers.altseg));
tcellids = cell(length(Markers.altseg));
trows = false(max(CellID),length(Markers.altseg));
%
for i1 = 1:length(Markers.altseg)
    markalt = Markers.altseg{i1};
    %
    % get folder and image names for altseg
    %
    fdname = [extractBefore(p.fname.folder,Markers.all{1}),...
        markalt,'\'];
    iname = [fdname,extractBefore(p.fname.name,'cell_seg_data.txt'),...
        'binary_seg_maps.tif'];
    %
    % get rows of altseg cells
    %
    trows(:,i1) = strcmp(p.fig.Phenotype, markalt);
    %
    % get inForm cellids of altseg cells
    %
    tcellids{i1} = double(p.fig.CellNum(trows(:,i1),:));
    %
    % read in the image for segmentation 
    %
    im = imread(iname,4);
    %
    % convert it to a linear index of labels
    %
    im = label2idx(im);
    im3{i1} = im(1,tcellids{i1});
end
%
%get filenames for 1ry seg images
%
iname = [p.fname.folder,'\',extractBefore(p.fname.name,...
    'cell_seg_data.txt'),'binary_seg_maps.tif'];
%
% get cellids of 1ry seg cells
%
trowsall = sum(trows,2) > 0;
cellids = double(p.fig.CellNum(~trowsall,:));
%
% read in image convert to linear index 
%
im2 = imread(iname,4);
im4 = label2idx(im2);
%
% get image dimensions
%
s = {size(im2)};
p.size.T = 2; p.size.B = s{1}(1)-1; p.size.L = 2; p.size.R = s{1}(2)-1;
%
% Remove non 1ry cells
%
im4 = im4(1,cellids);
%
% get expected size of cell vector
%
s2 = max(CellID);
%
% compile the segmentation images
%
im5 = cell(1,s2);
%
% first input all altsegs
%
for i1 = 1:length(Markers.altseg)
    im5(CellID(trows(:,1))) = im3{1};
end
%
% next input the 1ry segmentations
%
im5(CellID(~trowsall)) = im4;
%
% fill the cells and get them in the proper format
%
[obj,s] = OrganizeCells(p,s,s2,im5);
%
% remove others in altsegs
%
objt = obj(trowsall);
objt = cat(1,objt{:});
%
X = round(p.fig.CellXPos(orows));
Y = round(p.fig.CellYPos(orows));
Oth = sub2ind(s{1},Y,X);
rows = ismember(Oth,objt);
%
% count number of deleted others
%
p.flags.segclean = sum(rows == 1);
%
% remove those others from the output p.fig
%
p.fig(orows(rows),:) = [];
%
% remove those others from the CellID vector
%
CellID(orows(rows)) = [];
%
% remove those others from the cell objects variable and save to p.obj
%
p.obj = obj(CellID);
%
% add on those cells without membranes to the bottom of p.fig
%
nmr.CellNum = nmr.CellID; 
p.fig = vertcat(p.fig,nmr);
%
% get new and final CellIDs
%
CellID = (1:1:height(p.fig))';
p.fig.CellID = CellID;
%
end
%% function: OrganizeCells; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% close the cells on the edge of the tissue and organize the numbers
% inside the object cell to a circle so that the inside of the cell can
% be filled
%% --------------------------------------------------------------
%%
function [obj,s] = OrganizeCells(p,s,s2,im5)
%
% convert linear index to x,y
%
s = repmat(s,1,s2);
[y,x] = cellfun(@ind2sub,s,im5,'Uni',0);
%
% locate cells on edge of tissue & fill edges of cells on edge of image
%
X = cellfun(@min,x);
Y = cellfun(@min,y);
Y2 = cellfun(@max,y);
X2 = cellfun(@max,x);
rows = find(Y == p.size.T|Y2 == p.size.B|X == p.size.L|X2==p.size.R);
[x(rows),y(rows)] = cellfun(@(x,y)...
    filledges([x,y],p.size.T,p.size.B,p.size.L,p.size.R),...
    x(1,rows),y(1,rows),'Uni',0);
%
% find midpoints of the cell objects
%
mx = cellfun(@mean,x);
my = cellfun(@mean,y);
%
% get the angle of each point from the midpoint
%
ang = cellfun(@(x,y,mx,my)atan2(y-my,x-mx),...
    x,y,num2cell(mx),num2cell(my),'Uni',0);
%
% sorting the angles from the midpoint we can get indicies for a circle
%
[~,o] = cellfun(@sort,ang,'Uni',0);
%
% order the cells in the circle
%
x = cellfun(@(x,o)x(o),x,o,'Uni',0);
y = cellfun(@(y,o)y(o),y,o,'Uni',0);
%
% put the cells in reference coordinate rectangles
%
x1 = cellfun(@(x,y)x-y,x,num2cell(X),'Uni',0);
y1 = cellfun(@(x,y)x-y,y,num2cell(Y),'Uni',0);
%
% fill objects
%
obj = arrayfun(@(X,X2,Y,Y2)zeros(Y2-Y+1,X2-X+1),X,X2,Y,Y2,'Uni',0);
obj = cellfun(@fillobj,obj,y1,x1,'Uni',0);
[objy,objx] = cellfun(@(x)find(x==1),obj,'Uni',0);
%
% convert back to orginal image coordinates
%
objy = cellfun(@(x,y)x+y-1,objy,num2cell(Y),'Uni',0);
objx = cellfun(@(x,y)x+y-1,objx,num2cell(X),'Uni',0);
obj = cellfun(@sub2ind,s,objy,objx,'Uni',0);
%
% if the cell is very small then the obj is written as a column vector
% instead of a row vector
% to fix this find those cells
%
[~,r2] = cellfun(@size,obj,'Uni',0);
r2 = find(cat(1,r2{:})~=1);
%
% check if any of these cells exist
%
if ~isempty(r2)
    %
    % correct them to proper dimensions
    %
    for i4 = 1:length(r2)
        obj{r2(i4)} =  obj{r2(i4)}';
    end
end
end
%% function: filledges; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% this function takes in a cell object bon the edge of the image and 
% fills in the boundaries using T,B,L,R to define the edges as Top,
% Bottom, Left, Right
%% --------------------------------------------------------------
%%
function[x,y] = filledges(b,T,B,L,R)
w = 1;
tedge = [];
%
% get a vector that defines if the cell is on the top or bottom of the
% image
%
rows = b(:,2) == T |  b(:,2) == B;
edges = b(rows,:);
%
% find the edge where its pixel nieghborhood is missing a boundary
%
for i1 = 1:length(edges(:,1))
    x = edges(i1,:);
    neighbor=zeros(8,2);
    neighbor(1,:) = [x(1)-1, x(2)-1];
    neighbor(2,:) = [x(1),x(2)-1];
    neighbor(3,:) = [x(1)+1,x(2)-1];
    neighbor(4,:) = [x(1)-1,x(2)];
    neighbor(5,:) = [x(1)+1,x(2)];
    neighbor(6,:) = [x(1)-1,x(2)+1];
    neighbor(7,:) = [x(1),x(2)+1];
    neighbor(8,:) = [x(1)+1,x(2)+1];
    neighborhood(1) = sum(ismember(...
        neighbor([1,4,6],:),b,'rows'))>0;
    neighborhood(2) = sum(ismember(...
        neighbor([7,2],:),b,'rows'))>0;
    neighborhood(3) = sum(ismember(...
        neighbor([3,5,8],:),b,'rows'))>0;
    val = sum(neighborhood);
    if val == 1
        tedge(w,:) = x;
        w = w+1;
    end
end
%
% get a vector that defines if the cell is on the top or bottom of the
% image
%
rows = b(:,1) == L| b(:,1) ==  R;
edges = b(rows,:);
%
% find the edge where its pixel nieghborhood is missing a boundary
%
for i1 = 1:length(edges(:,1))
    x = edges(i1,:);
    neighbor=zeros(8,2);
    neighbor(1,:) = [x(1)-1, x(2)-1];
    neighbor(2,:) = [x(1),x(2)-1];
    neighbor(3,:) = [x(1)+1,x(2)-1];
    neighbor(4,:) = [x(1)-1,x(2)];
    neighbor(5,:) = [x(1)+1,x(2)];
    neighbor(6,:) = [x(1)-1,x(2)+1];
    neighbor(7,:) = [x(1),x(2)+1];
    neighbor(8,:) = [x(1)+1,x(2)+1];
    neighborhood(1) = sum(ismember(...
        neighbor([6,7,8],:),b,'rows'))>0;
    neighborhood(2) = sum(ismember(...
        neighbor([4,5],:),b,'rows'))>0;
    neighborhood(3) = sum(ismember(...
        neighbor([1,2,3],:),b,'rows'))>0;
    val = sum(neighborhood);
    if val == 1
        tedge(w,:) = x;
        w = w+1;
    end
end
%
% close the cell if it touches the edge
%
ns1 = [];ns2=[];
if ~isempty(tedge)&& length(tedge(:,1)) == 2 
    rows = tedge(:,1) == T | tedge(:,1) == B;
    %
    m1 = min(tedge(:,2))+1;
    m2 = max(tedge(:,2))-1;
    ns1(:,2) = m1:m2;
    %
    m1 = min(tedge(:,1))+1;
    m2 = max(tedge(:,1))-1;
    ns2(:,1) = m1:m2;
    %
    if sum(rows) == 0
        ns2(:,2) = tedge(1,2);
        ns = ns2;
    elseif sum(rows) == 2
        ns1(:,1) = tedge(1,1);
        ns = ns1;
    else
        ns1(:,1) = tedge(rows,1);
        ns2(:,2) = tedge(~rows,2);
        ns = vertcat(ns1,ns2);
    end
    %
    b = vertcat(b,ns);
    %
end
x = b(:,1);
y = b(:,2);
end
%% function: fillobj; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% fill the cell object 
%% --------------------------------------------------------------
%%
function[obj] = fillobj(obj,y,x)
s = size(obj);
x = x+1;
y = y+1;
obj(sub2ind(s,y,x)) = 1;
obj = imfill(obj);
end
%% function: getexprmark; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% determine the expression marker status of each cell; if the cell is in
% two cell membranes choose the cell that the expression marker call is
% closest to 
%% --------------------------------------------------------------
%%
function o = getexprmark(q,Markers)
o = q;
if ~isempty(o.fig)
    bin = [1,2,4,8,16,32,64,128];
    ExprPhenotype = zeros(height(o.fig),1);
    for i3 = 1:length(Markers.expr)
        mark = lower(Markers.expr{i3});
        % 
        % get cell indx
        %
        s = [o.size.B+1,o.size.R+1];
        ii = sub2ind(s,o.(mark).CellYPos,o.(mark).CellXPos);
        ii1 = zeros(1,height(o.fig))';
        xx1 = zeros(1,height(o.fig))';
        %
        [x2,x,rows] = findexpr(ii,o);
        %
        % assign cells which are inside two cells or outside of cell 
        % membranes; based on the cell centers distance
        %
        [~,x1] = min(edist(o.fig, o.(mark)(rows,:)),[],1);
        xx1(x1) = 1;
        %
        % for cells located only in one cell apply the expression to that cell
        %
        xx1(x2) = 1;
        %
        ExprPhenotype = ExprPhenotype + (xx1 .* bin(i3));
    end
    o.fig.ExprPhenotype = ExprPhenotype;
else
    o.fig.ExprPhenotype = zeros(0); 
end
%
end
%% function: findexpr; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% find expression marker calls which are inside two cells or outside of 
% cell membranes based on its cell center location and the mmembrane
% objects
%% --------------------------------------------------------------
%%
function[x2,x3,rows] = findexpr(ii,o)
%
% get the indicies of those rows in expression table that are inside at
% least one cell
%
x = cellfun(@(obj)find(ismember(ii,obj)),o.obj,'Uni',0);
x = cat(1,x{:});
% determine if that cell is in two cells or unique to just one
[a,b] = unique(x);
rows = ismember(x,x(setdiff(1:length(x),b)));

i2 = ii;
i2(a) = [];
i2 = [i2;ii(unique(x(rows)))];
% find those rows in expression table
rows = ismember(ii,i2);
%
% for cells located only in one cell apply the expression to that cell
%
ii = ii(~rows);
x3 = cellfun(@(obj)find(ismember(ii,obj)),o.obj,'Uni',0);
x2 = find(~cellfun(@isempty,x3))';
% if lineage cell has two cells keep the first one for now
x3 = cellfun(@(x)x(1),x3(x2));
end  
%% function: parsave; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% write the resulting table using writetable to the *\Table directory
% with the name *_cleaned_Phenotype_table.csv' and if the .csv file has
% more than 60 tumor cells and 600 cells total write the table for the
% ImageLoop for QC output
%% --------------------------------------------------------------
%%
function mktmp = parsave(fData,fname,Markers,wd)
%
% first check that there are enough columns to match the 9 color protocol
%
if width(fData.fig) < 83
    vars.Opals = {'DAPI','480','520','540','570','620','650','690','780'};
    vars.type = [repmat({'Mean'},1,4), repmat({'Total'},1,4)];
    vars.comp = repmat({'Nucleus', 'Membrane','EntireCell','Cytoplasm'},1,2);
    %
    for i3 = 1:length(vars.Opals)
        %
        b(:,i3) = cellfun(@(x,y)strcat(y,x,vars.Opals{i3}),...
            vars.comp,vars.type,'Uni',0);
        %
    end
    %    
    zz = b';
    zz = vertcat(zz(:))';
    vars.names = zz;
    %
    w = fData.fig;
    names_in = w.Properties.VariableNames;
    names_out = [{'CellID','SampleName','SlideID','fx','fy',...
        'CellNum','Phenotype','CellXPos',...
        'CellYPos','EntireCellArea'},vars.names{:},{'ExprPhenotype'}];
    ii = ~ismember(names_out,names_in);
    names_add = names_out(ii);
    vec_add = -1 * ones(height(w), 1); 
    %
    for i4 = 1:length(names_add)
        w.(names_add{i4}) = vec_add;
    end
    w = w(:,names_out);
    fData.fig = w;
elseif width(fData.fig) > 83
    disp('Warning: For database upload inForm output should contain ',...
        'no more than 9 color data');
end
%
% write out whole image csv table
%
nm = [wd,'\Results\Tables\',erase(fname.name,...
    'cell_seg_data.txt'),'cleaned_phenotype_table.csv'];
writetable(fData.fig(:,[1,3:end]),nm);
%
% write out image for use in figures later
% saves on images for figures if it is more than 600 cells and 60
% tumor cells in the image (if tumor marker exists)
%
if ~isempty(Markers.Tumor{1})
    r = strcmp(fData.fig.Phenotype,Markers.Tumor{1});
else
    r = ones(height(fData.fig));
end
if length(fData.fig.CellID) > 400 && length(find(r)) > 60
    fname2 = strcat(wd,'\Results\tmp_ForFiguresTables\',...
        erase(fname.name,'cell_seg_data.txt'),...
        'cleaned_phenotype_table.mat');
    save(fname2,'fData');
    mktmp = 'TRUE';
else
    mktmp = 'FALSE';
end
end
%% function: edist; 
%% --------------------------------------------------------------
%% Created by: Alex Szalay - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% compute the Euclidean distance between the two 2D vectors
%% --------------------------------------------------------------
%%
function e = edist(a,b)
    %
    e = zeros(height(a),height(b));
    for i=1:numel(a.CellXPos)
        X = double(a.CellXPos(i));
        Y = double(a.CellYPos(i));
        e(i,:) = sqrt(double(((b.CellXPos-X).^2+(b.CellYPos-Y).^2)));
    end
end
%%
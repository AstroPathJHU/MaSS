%% function: readtxt;
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description:
% reads in a single cell seg data file
%% --------------------------------------------------------
%%
function [Q, units] = readtxt(filnm,marker,wd,layers)
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
    { ' %s '},repmat('%f32 ',[1,4]),{' '},repmat('%f32 ',[1,5*layers]),...
    { ' %f32 %s '},repmat('%f32 ',[1,4]),{' '},repmat('%f32 ',[1,5*layers]),...
    { ' %f32 %s '},repmat('%f32 ',[1,4]),{' '},repmat('%f32 ',[1,5*layers]),...
     { ' %f32 %s '},repmat('%f32 ',[1,4]),{' '},repmat('%f32 ',[1,5*layers]),...
    {' '},repmat('%s ',[1,2]),{' '}, repmat('%f32 ',[1,4]),{' '}, ....
    repmat('%s ',[1,2]));
formatspec = formatspec{1};
%
try
    nm1 = [wd,'\Phenotyped\',marker,'\',filnm.name];
    T = readtable(nm1, 'Format', formatspec,...
        'Delimiter', '\t', 'TreatAsEmpty', {' ','#N/A'});
catch E
    disp(E);
    %
    nm = extractBefore(filnm.name,"]_cell_seg");
    if isempty(nm)
        nm = extractBefore(filnm.name,"]_CELL_SEG");
    end
    %
    try
        nm1 = [wd,'\Phenotyped\',marker,'\', nm, ']_CELL_SEG_DATA.TXT'];
        T = readtable(nm1,'Format',formatspec,...
            'Delimiter','\t','TreatAsEmpty',{' ','#N/A'});
    catch E
        disp(E);
         nm1 = [wd,'\Phenotyped\',marker,'\', nm, ']_cell_seg_data.txt'];
         T = readtable(nm1,'Format',formatspec,...
        'Delimiter','\t','TreatAsEmpty',{' ','#N/A'});
    end
    %
end
%
% Find Coordinates
%
f = strsplit(filnm.name,{'_[',']_'});
f = strsplit(f{2},',');
T.fx = repmat(str2double(f{1}),height(T),1);
T.fy = repmat(str2double(f{2}),height(T),1);
colnames = T.Properties.VariableNames;
%
if any(contains(colnames,'pixels'))
    units = 'pixels';
    areavariable = 'EntireCellArea_pixels_';
else
    units = 'microns';
    areavariable = 'EntireCellArea_squareMicrons_';
end
%
% to include autoflourescence add AutofluorescenceMean &
% AutofluorescenceTotal to the search and add ii{5} + ii{6} to 74
%
ii = cellfun(@(x) contains(colnames,x),...
    {'_Mean_', '_Total_', 'DAPIMean', 'DAPITotal'},...
    'UniformOutput', false);
ii3 = (ii{1} + ii{2} + ii{3} + ii{4}) > 0;
%
basecolumns = [{'SampleName','SlideID','fx','fy',...
    'CellID','Phenotype','CellXPosition',...
    'CellYPosition'}, areavariable, colnames(ii3)];
%
newnames = replace(basecolumns, '_NormalizedCounts_TotalWeighting_', '');
newnames = replace(newnames, 'Position', 'Pos');
OpalTypes = {'DAPI','480','520','540','570','620','650','690','780', 'AF'};
CellCompartments = {'Nucleus','Membrane','Cytoplasm','EntireCell'};
datatypes = {'Mean', 'Total'};
%
for i1 = 1:length(OpalTypes)
    for i2 = 1:length(CellCompartments)
        for i3 = 1:length(datatypes)
            ii = contains(newnames, OpalTypes{i1}) & ...
                contains(newnames, datatypes{i3}) & ...
                contains(newnames, CellCompartments{i2});
            if any(ii)
                newnames{ii} = [datatypes{i3},CellCompartments{i2},OpalTypes{i1}];
            end
        end
    end
end
%
Q = T(:, basecolumns);
Q.Properties.VariableNames = newnames;
%
end

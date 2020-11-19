function [Q, units] = readtxt(filnm,marker,vars,wd,layers)
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
try
    T = readtable([wd,'\Phenotyped\',marker,'\',filnm.name],'Format',formatspec,...
        'Delimiter','\t','TreatAsEmpty',{' ','#N/A'});
catch
    nm = extractBefore(filnm.name,"]_cell_seg");
    if isempty(nm)
        nm = extractBefore(filnm.name,"]_CELL_SEG");
    end
    try
        nm1 = [nm, ']_CELL_SEG_DATA.TXT'];
        T = readtable([wd,'\Phenotyped\',marker,'\',nm1],'Format',formatspec,...
            'Delimiter','\t','TreatAsEmpty',{' ','#N/A'});
    catch
         nm1 = [nm, ']_cell_seg_data.txt'];
         T = readtable([wd,'\Phenotyped\',marker,'\',nm1],'Format',formatspec,...
        'Delimiter','\t','TreatAsEmpty',{' ','#N/A'});
    end
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
% Selects variables of interest from originial data table by name
%
namestry(1,:) = [{'SampleName','SlideID','fx','fy',...
    'CellID','Phenotype','CellXPosition',...
    'CellYPosition'},areavariable,...
    vars.namesin{:}];
namestry(2,:) = [{'SampleName','SlideID','fx','fy',...
    'CellID','Phenotype','CellXPosition',...
    'CellYPosition'},areavariable,...
    vars.namesin_original{:}];
namestry(3,:) =  [{'SampleName','SlideID','fx','fy',...
    'CellID','Phenotype','CellXPosition',...
    'CellYPosition'},areavariable,...
    vars.namesin_dapierror{:}];
namestry(4,:) =  [{'SampleName','SlideID','fx','fy',...
    'CellID','Phenotype','CellXPosition',...
    'CellYPosition'},areavariable,...
    vars.namesin_dapierror_original{:}];
%
fail = 1;
icount = 1;
while fail == 1
    try
        Q = T(:,namestry(icount,:));
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
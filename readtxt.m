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

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
function mktmp = parsave(fData, fname, Markers, wd, imall)
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
    %
    % remove NaNs
    %
    for N = vars.names
        temp = num2cell(w.(N{1}));
        temp(cellfun(@isnan,temp)) = {[]};
        w.(N{1}) = temp;
    end
    fData.fig = w;
elseif width(fData.fig) > 83
    disp('Warning: For database upload inForm output should contain ',...
        'no more than 9 color data. Nulls not replaced in tables.');
end
%
% write out whole image csv table
%
nm = [wd,'\Phenotyped\Results\Tables\',extractBefore(fname.name,...
    ']_cell_seg'),']_cleaned_phenotype_table.csv'];
if isempty(nm)
    nm = [wd,'\Phenotyped\Results\Tables\',extractBefore(fname.name,...
    ']_CELL_SEG'),']_cleaned_phenotype_table.csv'];
end
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
if (length(fData.fig.CellID) > 400) && (imall || (length(find(r)) > 60))
    fname2 = strcat(wd,'\Phenotyped\Results\tmp_ForFiguresTables\',...
        extractBefore(fname.name,']_cell_seg'),...
        ']_cleaned_phenotype_table.mat');
    if isempty(fname2)
        fname2 = strcat(wd,'\Phenotyped\Results\tmp_ForFiguresTables\',...
            extractBefore(fname.name,']_CELL_SEG'),...
            ']_cleaned_phenotype_table.mat');
    end
    save(fname2,'fData');
    mktmp = 'TRUE';
else
    mktmp = 'FALSE';
end
end

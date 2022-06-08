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

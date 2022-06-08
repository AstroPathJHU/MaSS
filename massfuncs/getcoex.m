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

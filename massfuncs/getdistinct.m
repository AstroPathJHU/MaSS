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

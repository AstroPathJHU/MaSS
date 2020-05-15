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
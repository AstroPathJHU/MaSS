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
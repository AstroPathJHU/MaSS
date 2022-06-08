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

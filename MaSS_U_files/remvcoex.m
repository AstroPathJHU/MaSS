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
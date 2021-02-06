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
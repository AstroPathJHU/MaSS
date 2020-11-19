function o = getexprmark(q,Markers)
o = q;
if ~isempty(o.fig)
    ii = ismember([480;520;540;570;620;650;690;780],Markers.Opals);
    bins = [2,4,8,16,32,64,128,256];
    bins(~ii) = [];
    ii = ismember(Markers.all,Markers.expr);
    bins = bins(ii);
    %
    ExprPhenotype = zeros(height(o.fig),1);
    for i3 = 1:length(Markers.expr)
        mark = lower(Markers.expr{i3});
        % 
        % get cell indx
        %
        s = [o.size.B+1,o.size.R+1];
        mx = o.(mark).CellXPos + 1;
        my = o.(mark).CellYPos + 1;
        ii = sub2ind(s, my, mx);
        xx1 = zeros(1,height(o.fig))';
        %
        [x2,x,rows] = findexpr(ii,o);
        %
        % assign cells which are inside two cells or outside of cell 
        % membranes; based on the cell centers distance
        %
        [~,x1] = min(edist(o.fig, o.(mark)(rows,:)),[],1);
        xx1(x1) = 1;
        %
        % for cells located only in one cell apply the expression to that cell
        %
        xx1(x2) = 1;
        %
        ExprPhenotype = ExprPhenotype + (xx1 .* bins(i3));
    end
    o.fig.ExprPhenotype = ExprPhenotype;
else
    o.fig.ExprPhenotype = zeros(0); 
end
%
end

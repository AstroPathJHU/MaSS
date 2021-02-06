function e = edist(a,b)
    %
    e = zeros(height(a),height(b));
    for i=1:numel(a.CellXPos)
        X = double(a.CellXPos(i));
        Y = double(a.CellYPos(i));
        e(i,:) = sqrt(double(((b.CellXPos-X).^2+(b.CellYPos-Y).^2)));
    end
end

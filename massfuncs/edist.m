%% function: edist; 
%% --------------------------------------------------------------
%% Created by: Alex Szalay - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% compute the Euclidean distance between the two 2D vectors
%% --------------------------------------------------------------
%%
function e = edist(a,b)
    %
    e = zeros(height(a),height(b));
    for i=1:numel(a.CellXPos)
        X = double(a.CellXPos(i));
        Y = double(a.CellYPos(i));
        e(i,:) = sqrt(double(((b.CellXPos-X).^2+(b.CellYPos-Y).^2)));
    end
end

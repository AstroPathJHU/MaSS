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

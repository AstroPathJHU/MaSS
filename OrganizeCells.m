function [obj,s] = OrganizeCells(p,s,s2,im5)
%
% convert linear index to x,y
%
s = repmat(s,1,s2);
[y,x] = cellfun(@ind2sub,s,im5,'Uni',0);
%
% locate cells on edge of tissue & fill edges of cells on edge of image
%
X = cellfun(@min,x);
Y = cellfun(@min,y);
Y2 = cellfun(@max,y);
X2 = cellfun(@max,x);
rows = find(Y == p.size.T|Y2 == p.size.B|X == p.size.L|X2==p.size.R);
[x(rows),y(rows)] = cellfun(@(x,y)...
    filledges([x,y],p.size.T,p.size.B,p.size.L,p.size.R),...
    x(1,rows),y(1,rows),'Uni',0);
%
% find midpoints of the cell objects
%
mx = cellfun(@mean,x);
my = cellfun(@mean,y);
%
% get the angle of each point from the midpoint
%
ang = cellfun(@(x,y,mx,my)atan2(y-my,x-mx),...
    x,y,num2cell(mx),num2cell(my),'Uni',0);
%
% sorting the angles from the midpoint we can get indicies for a circle
%
[~,o] = cellfun(@sort,ang,'Uni',0);
%
% order the cells in the circle
%
x = cellfun(@(x,o)x(o),x,o,'Uni',0);
y = cellfun(@(y,o)y(o),y,o,'Uni',0);
%
% put the cells in reference coordinate rectangles
%
x1 = cellfun(@(x,y)x-y,x,num2cell(X),'Uni',0);
y1 = cellfun(@(x,y)x-y,y,num2cell(Y),'Uni',0);
%
% fill objects
%
obj = arrayfun(@(X,X2,Y,Y2)zeros(Y2-Y+1,X2-X+1),X,X2,Y,Y2,'Uni',0);
obj = cellfun(@fillobj,obj,y1,x1,'Uni',0);
[objy,objx] = cellfun(@(x)find(x==1),obj,'Uni',0);
%
% convert back to orginal image coordinates
%
objy = cellfun(@(x,y)x+y-1,objy,num2cell(Y),'Uni',0);
objx = cellfun(@(x,y)x+y-1,objx,num2cell(X),'Uni',0);
obj = cellfun(@sub2ind,s,objy,objx,'Uni',0);
%
% if the cell is very small then the obj is written as a column vector
% instead of a row vector
% to fix this find those cells
%
[~,r2] = cellfun(@size,obj,'Uni',0);
r2 = find(cat(1,r2{:})~=1);
%
% check if any of these cells exist
%
if ~isempty(r2)
    %
    % correct them to proper dimensions
    %
    for i4 = 1:length(r2)
        obj{r2(i4)} =  obj{r2(i4)}';
    end
end
end
%% function: filledges; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% this function takes in a cell object bon the edge of the image and 
% fills in the boundaries using T,B,L,R to define the edges as Top,
% Bottom, Left, Right
%% --------------------------------------------------------------
%%
function[x,y] = filledges(b,T,B,L,R)
w = 1;
tedge = [];
%
% get a vector that defines if the cell is on the top or bottom of the
% image
%
rows = b(:,2) == T |  b(:,2) == B;
edges = b(rows,:);
%
% find the edge where its pixel nieghborhood is missing a boundary
%
for i1 = 1:length(edges(:,1))
    x = edges(i1,:);
    neighbor=zeros(8,2);
    neighbor(1,:) = [x(1)-1, x(2)-1];
    neighbor(2,:) = [x(1),x(2)-1];
    neighbor(3,:) = [x(1)+1,x(2)-1];
    neighbor(4,:) = [x(1)-1,x(2)];
    neighbor(5,:) = [x(1)+1,x(2)];
    neighbor(6,:) = [x(1)-1,x(2)+1];
    neighbor(7,:) = [x(1),x(2)+1];
    neighbor(8,:) = [x(1)+1,x(2)+1];
    neighborhood(1) = sum(ismember(...
        neighbor([1,4,6],:),b,'rows'))>0;
    neighborhood(2) = sum(ismember(...
        neighbor([7,2],:),b,'rows'))>0;
    neighborhood(3) = sum(ismember(...
        neighbor([3,5,8],:),b,'rows'))>0;
    val = sum(neighborhood);
    if val == 1
        tedge(w,:) = x;
        w = w+1;
    end
end
%
% get a vector that defines if the cell is on the top or bottom of the
% image
%
rows = b(:,1) == L| b(:,1) ==  R;
edges = b(rows,:);
%
% find the edge where its pixel nieghborhood is missing a boundary
%
for i1 = 1:length(edges(:,1))
    x = edges(i1,:);
    neighbor=zeros(8,2);
    neighbor(1,:) = [x(1)-1, x(2)-1];
    neighbor(2,:) = [x(1),x(2)-1];
    neighbor(3,:) = [x(1)+1,x(2)-1];
    neighbor(4,:) = [x(1)-1,x(2)];
    neighbor(5,:) = [x(1)+1,x(2)];
    neighbor(6,:) = [x(1)-1,x(2)+1];
    neighbor(7,:) = [x(1),x(2)+1];
    neighbor(8,:) = [x(1)+1,x(2)+1];
    neighborhood(1) = sum(ismember(...
        neighbor([6,7,8],:),b,'rows'))>0;
    neighborhood(2) = sum(ismember(...
        neighbor([4,5],:),b,'rows'))>0;
    neighborhood(3) = sum(ismember(...
        neighbor([1,2,3],:),b,'rows'))>0;
    val = sum(neighborhood);
    if val == 1
        tedge(w,:) = x;
        w = w+1;
    end
end
%
% close the cell if it touches the edge
%
ns1 = [];ns2=[];
if ~isempty(tedge)&& length(tedge(:,1)) == 2 
    rows = tedge(:,1) == T | tedge(:,1) == B;
    %
    m1 = min(tedge(:,2))+1;
    m2 = max(tedge(:,2))-1;
    ns1(:,2) = m1:m2;
    %
    m1 = min(tedge(:,1))+1;
    m2 = max(tedge(:,1))-1;
    ns2(:,1) = m1:m2;
    %
    if sum(rows) == 0
        ns2(:,2) = tedge(1,2);
        ns = ns2;
    elseif sum(rows) == 2
        ns1(:,1) = tedge(1,1);
        ns = ns1;
    else
        ns1(:,1) = tedge(rows,1);
        ns2(:,2) = tedge(~rows,2);
        ns = vertcat(ns1,ns2);
    end
    %
    b = vertcat(b,ns);
    %
end
x = b(:,1);
y = b(:,2);
end
%% function: fillobj; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% fill the cell object 
%% --------------------------------------------------------------
%%
function[obj] = fillobj(obj,y,x)
s = size(obj);
x = x+1;
y = y+1;
obj(sub2ind(s,y,x)) = 1;
obj = imfill(obj);
end
%% function: getexprmark; 
%% --------------------------------------------------------------
%% Created by: Benjamin Green - Johns Hopkins - 01/03/2019
%% --------------------------------------------------------------
%% Description
% determine the expression marker status of each cell; if the cell is in
% two cell membranes choose the cell that the expression marker call is
% closest to 
%% --------------------------------------------------------------
%%
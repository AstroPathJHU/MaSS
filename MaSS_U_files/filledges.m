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
